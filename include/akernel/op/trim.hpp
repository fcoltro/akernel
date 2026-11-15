// akernel — trimmed curves, trimmed-face loops and certified point
// classification.
//
// A trimmed face is a carrier patch plus closed loops in its parameter
// square. A loop is made of straight pieces along the sides of the square
// and of arcs, which are parts of certified intersection curves
// (solve/ssi.hpp) seen in this face's two parameters.
//
// Classification of a parameter point q against a loop is by winding number,
// and it is exact, not approximate:
//   - every piece of the true loop is enclosed in a convex box together with
//     a straight chord between two points of that box;
//   - if q lies outside the box, the true piece and its chord are homotopic
//     in the plane minus q, so replacing the piece by the chord does not
//     change the winding number;
//   - the winding number of the resulting polygon is computed with exact
//     orientation predicates.
// If q lies inside a box, the box is refined (the arc enclosure is
// recomputed for half the slice range). If that does not separate q from the
// curve within the depth limit, the answer is "undecided", never a guess.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

#include <akernel/num/predicates.hpp>
#include <akernel/solve/ssi.hpp>

namespace ak::op {

using num::Interval;
using solve::Box4;
using solve::Pt4;
using solve::SsiCurve;
using solve::SsiSystem;
using solve::TubeSegment;

/// Slack allowed between the nominal slice value of a join and the slice of
/// the exact point there. Joins are verified against it before any
/// refinement is allowed.
inline constexpr double kSlicePad = 1e-9;

enum class Tri { Yes, No, Unknown };

/// Affine map between two parameter planes, with interval coefficients:
///   u' = a[0][0] + a[0][1] u + a[0][2] v,   v' = a[1][0] + a[1][1] u + a[1][2] v.
/// Used between two planar faces that lie in the same plane; the coefficients
/// enclose the exact map, so mapped boxes enclose the exact images.
/// An exactly representable value of an expansion, if it has one.
inline bool exact_value(const num::Expansion& e, double& x) {
    if (e.size() > 1) return false;
    x = e.estimate();
    return true;
}
inline bool is_point(const Interval& x) { return x.lo == x.hi; }

/// c0 + c1 x + c2 y. When every input is a point and the exact result is a
/// double, the result is that point: interval arithmetic would widen it by a
/// unit even when nothing was rounded, and an exact parameter (a vertex at a
/// corner, say) must stay exact through a dyadic map for the exact tests that
/// use it.
inline Interval affine_combination(const Interval& c0, const Interval& c1, const Interval& x, const Interval& c2,
                                   const Interval& y) {
    if (is_point(c0) && is_point(c1) && is_point(x) && is_point(c2) && is_point(y)) {
        using E = num::Expansion;
        double r = 0.0;
        if (exact_value(E(c0.lo) + E(c1.lo) * E(x.lo) + E(c2.lo) * E(y.lo), r)) return Interval(r);
    }
    return c0 + c1 * x + c2 * y;
}

struct Affine {
    Interval a[2][3];
    static Affine identity() {
        Affine m;
        m.a[0][0] = Interval(0.0); m.a[0][1] = Interval(1.0); m.a[0][2] = Interval(0.0);
        m.a[1][0] = Interval(0.0); m.a[1][1] = Interval(0.0); m.a[1][2] = Interval(1.0);
        return m;
    }
    void apply(const Interval& u, const Interval& v, Interval& uo, Interval& vo) const {
        uo = affine_combination(a[0][0], a[0][1], u, a[0][2], v);
        vo = affine_combination(a[1][0], a[1][1], u, a[1][2], v);
    }
    /// The linear part, for directions.
    void apply_linear(const Interval& du, const Interval& dv, Interval& uo, Interval& vo) const {
        const Interval zero(0.0);
        uo = affine_combination(zero, a[0][1], du, a[0][2], dv);
        vo = affine_combination(zero, a[1][1], du, a[1][2], dv);
    }
    Interval det() const { return a[0][1] * a[1][2] - a[0][2] * a[1][1]; }
};

// Travel direction along a curve is from its front to its back. A segment
// with `reversed` set is travelled from b to a.
inline double seg_start_slice(const TubeSegment& s) { return s.reversed ? s.len : 0.0; }
inline double seg_end_slice(const TubeSegment& s) { return s.reversed ? 0.0 : s.len; }
inline const Pt4& seg_start_point(const TubeSegment& s) { return s.reversed ? s.b : s.a; }
inline const Pt4& seg_end_point(const TubeSegment& s) { return s.reversed ? s.a : s.b; }

/// Rigorous slice value of a box in a segment's frame.
inline Interval slice_of(const Box4& P, const TubeSegment& s) {
    Interval c(0.0);
    for (int i = 0; i < 4; ++i) c += Interval(s.m[i]) * (P[i] - Interval(s.a[i]));
    return c;
}

/// Is the box proven inside (Yes), proven outside (No) the certified region
/// of a segment, or neither?
inline Tri region_status(const Box4& P, const TubeSegment& s) {
    if (solve::in_region(P, s)) return Tri::Yes;
    for (int i = 0; i < 4; ++i)
        if (num::disjoint(P[i], s.X[i])) return Tri::No;
    const Interval c = slice_of(P, s);
    if (c.hi < s.c_lo || c.lo > s.c_hi) return Tri::No;
    return Tri::Unknown;
}

struct Span {
    int seg = 0;
    double ca = 0.0, cb = 0.0; ///< slice range in travel order
};

/// A whole certified intersection curve of two carriers, with its travel
/// spans. A closed curve has one more span than segments: the closing span
/// on segment 0, from where the last segment ends to where the first starts.
struct CurveBase {
    std::shared_ptr<const SsiSystem> sys;
    std::shared_ptr<const SsiCurve> curve;
    bool closed = false;
    std::vector<Span> spans;
    std::vector<char> join_ok; ///< join between span j and j+1 verified to within kSlicePad
};

inline std::shared_ptr<CurveBase> make_curve_base(std::shared_ptr<const SsiSystem> sys,
                                                  std::shared_ptr<const SsiCurve> curve) {
    auto base = std::make_shared<CurveBase>();
    base->sys = std::move(sys);
    base->curve = std::move(curve);
    base->closed = base->curve->closed;
    const auto& segs = base->curve->segs;
    const int n = static_cast<int>(segs.size());
    for (int k = 0; k < n; ++k) base->spans.push_back({k, seg_start_slice(segs[k]), seg_end_slice(segs[k])});
    // The exact end point of `prev` must sit within the pad of the nominal
    // slice at which `next` takes over.
    auto join = [&](const TubeSegment& prev, const TubeSegment& next, double next_slice) {
        Pt4 x = seg_end_point(prev);
        Box4 enc;
        if (!solve::certify_point(*base->sys, prev.m, prev.a, seg_end_slice(prev), x, enc)) return false;
        const Interval c = slice_of(enc, next);
        return next_slice - kSlicePad <= c.lo && c.hi <= next_slice + kSlicePad;
    };
    for (int k = 1; k < n; ++k) base->join_ok.push_back(join(segs[k - 1], segs[k], seg_start_slice(segs[k])) ? 1 : 0);
    if (base->closed) {
        const TubeSegment& first = segs.front();
        const Pt4& e = seg_end_point(segs.back());
        double ce = 0.0;
        for (int i = 0; i < 4; ++i) ce += first.m[i] * (e[i] - first.a[i]);
        base->spans.push_back({0, ce, seg_start_slice(first)});
        base->join_ok.push_back(join(segs.back(), first, ce) ? 1 : 0);
    }
    return base;
}

/// A certified position on a curve: a span and a slice interval in the frame
/// of that span's segment.
struct CurvePos {
    int span = 0;
    Interval c;
    bool dec = false; ///< the span is travelled with decreasing slice value
    Interval travel() const { return dec ? -c : c; }
};

inline CurvePos make_pos(const CurveBase& B, int span, const Interval& c) {
    const Span& s = B.spans[static_cast<std::size_t>(span)];
    return {span, c, s.cb < s.ca};
}

/// Is p strictly before q in the travel order of the span list?
inline Tri before(const CurvePos& p, const CurvePos& q) {
    if (p.span != q.span) return p.span < q.span ? Tri::Yes : Tri::No;
    const Interval a = p.travel(), b = q.travel();
    if (a.hi < b.lo) return Tri::Yes;
    if (b.hi < a.lo) return Tri::No;
    return Tri::Unknown;
}

/// Locate a certified zero of the curve's system (given by its enclosure P)
/// on the curve. Returns false if that could not be decided. `on` says
/// whether the point lies on the curve at all.
///
/// A position is only ever given in a span whose slice range it meets. A
/// segment's certified region also covers margins beyond that range (its
/// overlap with the neighbouring segments, and beyond the two ends of an
/// open curve). A point in an inner margin lies in the neighbouring span's
/// range and is found there; a point that lies wholly in a margin beyond an
/// end of the curve is on the continuation of the polynomial curve outside
/// the parameter domain, not on the curve. A point whose enclosure meets the
/// range only at its end (a vertex where the curve leaves a face, known to
/// rounding) is placed in that span.
inline bool locate_on_curve(const CurveBase& B, const Box4& P, bool& on, CurvePos& pos) {
    const auto& segs = B.curve->segs;
    on = false;
    std::vector<signed char> status(segs.size(), -1); // -1 unknown yet, 0 outside, 1 inside
    for (int j = 0; j < static_cast<int>(B.spans.size()); ++j) {
        const Span& sp = B.spans[static_cast<std::size_t>(j)];
        signed char& st = status[static_cast<std::size_t>(sp.seg)];
        if (st < 0) {
            const Tri t = region_status(P, segs[static_cast<std::size_t>(sp.seg)]);
            if (t == Tri::Unknown) return false;
            st = t == Tri::Yes ? 1 : 0;
        }
        if (st == 0) continue;
        const CurvePos p = make_pos(B, j, slice_of(P, segs[static_cast<std::size_t>(sp.seg)]));
        const Interval t = p.travel();
        double s0 = p.dec ? -sp.ca : sp.ca, s1 = p.dec ? -sp.cb : sp.cb;
        // At the two ends of an open curve, a vertex where the curve leaves a
        // parameter square can come out a few units beyond the end in slice
        // value (the slice is measured with the segment's rounded direction).
        // Allow the join slack there; a point beyond the end by more is on
        // the continuation outside the square, not on the curve.
        if (!B.closed && j == 0) s0 -= kSlicePad;
        if (!B.closed && j + 1 == static_cast<int>(B.spans.size())) s1 += kSlicePad;
        if (t.hi < s0 || t.lo > s1) continue; // wholly in a margin of this segment
        on = true;
        pos = p;
        return true;
    }
    return true;
}

/// Part of a curve between two certified positions: the geometry of an edge.
struct ArcGeom {
    std::shared_ptr<const CurveBase> base;
    CurvePos a, b;     ///< ends, in travel order
    bool wrap = false; ///< runs past the end of the span list and on from its start
    std::vector<Span> spans;
    bool refinable = false; ///< ends and inner joins verified to within kSlicePad
};

inline std::shared_ptr<ArcGeom> make_arc(std::shared_ptr<const CurveBase> base, const CurvePos& a,
                                         const CurvePos& b, bool wrap) {
    auto arc = std::make_shared<ArcGeom>();
    arc->base = std::move(base);
    arc->a = a;
    arc->b = b;
    arc->wrap = wrap;
    const CurveBase& B = *arc->base;
    const int n = static_cast<int>(B.spans.size());
    bool joins = true;
    auto piece = [&](int j, bool from_a, bool to_b) {
        Span s = B.spans[static_cast<std::size_t>(j)];
        if (from_a) s.ca = num::mid(a.c);
        if (to_b) s.cb = num::mid(b.c);
        arc->spans.push_back(s);
    };
    if (!wrap && a.span == b.span) {
        piece(a.span, true, true);
    } else {
        int j = a.span;
        piece(j, true, false);
        for (int guard = 0; guard <= n; ++guard) {
            // The join from the last span back to the first is the same
            // point of the same segment, so it needs no check.
            if (j + 1 < n && !B.join_ok[static_cast<std::size_t>(j)]) joins = false;
            j = (j + 1) % n;
            if (j == b.span) { piece(j, false, true); break; }
            piece(j, false, false);
        }
    }
    arc->refinable = joins && num::width(a.c) <= kSlicePad && num::width(b.c) <= kSlicePad;
    return arc;
}

/// The whole of a closed curve, from its seam to its seam.
inline std::shared_ptr<ArcGeom> make_full_arc(std::shared_ptr<const CurveBase> base) {
    const double s0 = base->spans.front().ca;
    const CurvePos a = make_pos(*base, 0, Interval(s0));
    const CurvePos b = make_pos(*base, static_cast<int>(base->spans.size()) - 1, Interval(s0));
    return make_arc(std::move(base), a, b, false);
}

/// A certified point of the arc strictly between its ends.
inline bool arc_point(const ArcGeom& arc, Box4& enc) {
    const Span& sp = arc.spans[arc.spans.size() / 2];
    const TubeSegment& seg = arc.base->curve->segs[static_cast<std::size_t>(sp.seg)];
    const double c = 0.5 * (sp.ca + sp.cb);
    Pt4 x;
    for (int i = 0; i < 4; ++i) x[i] = seg.a[i] + (c / seg.len) * (seg.b[i] - seg.a[i]);
    return solve::certify_point(*arc.base->sys, seg.m, seg.a, c, x, enc) && solve::in_region(enc, seg);
}

/// Rigorous enclosure of the arc of a segment over the slice values
/// [c0, c1], which must lie inside the segment's certified range. The
/// Krawczyk operator is used as a contractor: it encloses every zero in the
/// box it is given, and the segment is known to hold exactly one per slice.
inline Box4 arc_enclosure(const SsiSystem& sys, const TubeSegment& seg, double c0, double c1) {
    Box4 X = seg.X;
    const double cm = 0.5 * (c0 + c1);
    Pt4 xhat;
    for (int i = 0; i < 4; ++i) xhat[i] = seg.a[i] + (cm / seg.len) * (seg.b[i] - seg.a[i]);
    {
        Pt4 y = xhat;
        bool ok = solve::newton_slice(sys, seg.m, seg.a, cm, y);
        for (int i = 0; i < 4 && ok; ++i) ok = num::contains(X[i], y[i]);
        if (ok) xhat = y;
    }
    const Interval c(c0, c1);
    for (int it = 0; it < 4; ++it) {
        for (int i = 0; i < 4; ++i) xhat[i] = std::min(std::max(xhat[i], X[i].lo), X[i].hi);
        Box4 K;
        if (!solve::krawczyk_operator(sys, seg.m, seg.a, c, xhat, X, K)) break;
        Box4 Y;
        bool empty = false;
        for (int i = 0; i < 4; ++i) {
            const double lo = std::max(X[i].lo, K[i].lo), hi = std::min(X[i].hi, K[i].hi);
            if (!(lo <= hi)) empty = true;
            Y[i] = Interval(lo, hi);
        }
        if (empty) break; // cannot happen for a certified segment; keep the last enclosure
        // A fixed point: every further iteration would repeat this one exactly.
        bool same = true;
        for (int i = 0; i < 4 && same; ++i) same = Y[i].lo == X[i].lo && Y[i].hi == X[i].hi;
        if (same) break;
        X = Y;
    }
    return X;
}

/// Enclosure of the arc of a segment over a slice range, clamped to the
/// segment's certified range and padded so that the exact join points are
/// covered.
inline Box4 padded_enclosure(const SsiSystem& sys, const TubeSegment& seg, double ca, double cb) {
    const double lo = std::max(std::min(ca, cb) - kSlicePad, seg.c_lo);
    const double hi = std::min(std::max(ca, cb) + kSlicePad, seg.c_hi);
    return lo <= hi ? arc_enclosure(sys, seg, lo, hi) : seg.X;
}

/// A piece of a loop: a chord and a box that contains both the chord and the
/// true piece it stands for.
struct Link {
    double p0[2] = {0.0, 0.0}, p1[2] = {0.0, 0.0};
    Interval bu, bv;
    const ArcGeom* arc = nullptr; ///< null for a straight piece
    int coord = 0;                ///< 0: first carrier of the curve (u,v); 2: second (s,t)
    int seg = -1;
    double ca = 0.0, cb = 0.0;    ///< slice range in travel order
    int depth = 0;
    bool mapped = false;          ///< the curve's coordinates are mapped by `map` into this face's
    Affine map;
    bool exact_ends = false;      ///< straight piece whose end points are exact: the chord is the piece
    Interval e0[2], e1[2];        ///< straight piece: boxes holding its true end points and the chord's
};

/// A curve point in this face's coordinates.
inline void to_face(const Link& L, const Interval& u, const Interval& v, Interval& uo, Interval& vo) {
    if (L.mapped) L.map.apply(u, v, uo, vo);
    else { uo = u; vo = v; }
}

struct CycleGeom {
    std::vector<Link> links;
};

inline Interval hull_with(const Interval& b, double x, double y) {
    return Interval(std::min(b.lo, std::min(x, y)), std::max(b.hi, std::max(x, y)));
}

/// Straight piece between two nodes, each given with its enclosure.
inline void append_straight_link(std::vector<Link>& out, const double p0[2], const Interval box0[2],
                                 const double p1[2], const Interval box1[2]) {
    Link L;
    L.p0[0] = p0[0]; L.p0[1] = p0[1];
    L.p1[0] = p1[0]; L.p1[1] = p1[1];
    L.bu = hull_with(num::hull(box0[0], box1[0]), p0[0], p1[0]);
    L.bv = hull_with(num::hull(box0[1], box1[1]), p0[1], p1[1]);
    L.exact_ends = box0[0].lo == box0[0].hi && box0[1].lo == box0[1].hi && box1[0].lo == box1[0].hi &&
                   box1[1].lo == box1[1].hi;
    for (int i = 0; i < 2; ++i) {
        L.e0[i] = hull_with(box0[i], p0[i], p0[i]);
        L.e1[i] = hull_with(box1[i], p1[i], p1[i]);
    }
    out.push_back(L);
}

/// Is the query box clear of a straight piece whose ends are enclosures? The
/// true piece and the chord are both segments from a point of one end box to
/// a point of the other, so both lie in the convex hull of the two boxes, and
/// the chord can be moved onto the piece inside that hull. A strict sign of
/// (b - a) x (q - a) over every a, b in the end boxes and q in the query box
/// proves that no such segment passes through the query box.
inline bool clear_of_straight(const Link& L, const Interval& qu, const Interval& qv) {
    const Interval c = (L.e1[0] - L.e0[0]) * (qv - L.e0[1]) - (L.e1[1] - L.e0[1]) * (qu - L.e0[0]);
    return c.lo > 0.0 || c.hi < 0.0;
}

/// Links of an arc, travelled forward (a to b) or backward. pf and pb are
/// the chord end points to use at its two ends, in this face's coordinates.
/// With `map`, the curve's coordinates `coord` are mapped into the face's.
inline void append_arc_links(std::vector<Link>& out, const ArcGeom& arc, int coord, bool forward,
                             const double pf[2], const double pb[2], const Affine* map = nullptr) {
    const int ns = static_cast<int>(arc.spans.size());
    const std::size_t co = static_cast<std::size_t>(coord);
    const auto& segs = arc.base->curve->segs;
    Link proto;
    if (map) { proto.mapped = true; proto.map = *map; }
    // Chord points: P[0] = front, P[j] = join after span j-1, P[ns] = back.
    std::vector<std::array<double, 2>> P(static_cast<std::size_t>(ns + 1));
    P.front() = {pf[0], pf[1]};
    P.back() = {pb[0], pb[1]};
    for (int j = 1; j < ns; ++j) {
        const Span& sp = arc.spans[static_cast<std::size_t>(j - 1)];
        const TubeSegment& s = segs[static_cast<std::size_t>(sp.seg)];
        const Pt4& x = sp.cb == 0.0 ? s.a : s.b;
        Interval u, v;
        to_face(proto, Interval(x[co]), Interval(x[co + 1]), u, v);
        P[static_cast<std::size_t>(j)] = {num::mid(u), num::mid(v)};
    }
    for (int n = 0; n < ns; ++n) {
        const int j = forward ? n : ns - 1 - n;
        const Span& sp = arc.spans[static_cast<std::size_t>(j)];
        const TubeSegment& s = segs[static_cast<std::size_t>(sp.seg)];
        const auto& a = P[static_cast<std::size_t>(forward ? j : j + 1)];
        const auto& b = P[static_cast<std::size_t>(forward ? j + 1 : j)];
        Link L = proto;
        L.p0[0] = a[0]; L.p0[1] = a[1];
        L.p1[0] = b[0]; L.p1[1] = b[1];
        Interval bu, bv;
        to_face(proto, s.X[co], s.X[co + 1], bu, bv);
        L.bu = hull_with(bu, a[0], b[0]);
        L.bv = hull_with(bv, a[1], b[1]);
        L.arc = &arc;
        L.coord = coord;
        L.seg = sp.seg;
        L.ca = forward ? sp.ca : sp.cb;
        L.cb = forward ? sp.cb : sp.ca;
        out.push_back(L);
    }
}

/// Recompute a link's box from the arc enclosure of its slice range.
inline void tighten_link(Link& L) {
    const TubeSegment& seg = L.arc->base->curve->segs[static_cast<std::size_t>(L.seg)];
    const Box4 X = padded_enclosure(*L.arc->base->sys, seg, L.ca, L.cb);
    const std::size_t co = static_cast<std::size_t>(L.coord);
    Interval bu, bv;
    to_face(L, X[co], X[co + 1], bu, bv);
    L.bu = hull_with(bu, L.p0[0], L.p1[0]);
    L.bv = hull_with(bv, L.p0[1], L.p1[1]);
}

/// Split a link at the middle of its slice range.
inline void split_link(const Link& L, Link& first, Link& second) {
    const TubeSegment& seg = L.arc->base->curve->segs[static_cast<std::size_t>(L.seg)];
    const double cm = 0.5 * (L.ca + L.cb);
    Pt4 x;
    for (int i = 0; i < 4; ++i) x[i] = seg.a[i] + (cm / seg.len) * (seg.b[i] - seg.a[i]);
    {
        Pt4 y = x;
        bool ok = solve::newton_slice(*L.arc->base->sys, seg.m, seg.a, cm, y);
        for (int i = 0; i < 4 && ok; ++i) ok = std::isfinite(y[i]);
        if (ok) x = y;
    }
    const std::size_t co = static_cast<std::size_t>(L.coord);
    Interval mu, mv;
    to_face(L, Interval(x[co]), Interval(x[co + 1]), mu, mv);
    first = L;
    second = L;
    first.cb = cm;
    first.p1[0] = num::mid(mu); first.p1[1] = num::mid(mv);
    second.ca = cm;
    second.p0[0] = num::mid(mu); second.p0[1] = num::mid(mv);
    first.depth = second.depth = L.depth + 1;
    tighten_link(first);
    tighten_link(second);
}

struct WindingResult {
    bool decided = false;
    bool on_boundary = false; ///< undecided because the point is on a straight piece
    int wn = 0;
};

/// Winding number of the true loop about every point of the query box
/// (a point query is a degenerate box). Refines the loop's links as needed.
inline WindingResult winding(CycleGeom& cyc, const Interval& qu, const Interval& qv, int max_depth) {
    WindingResult res;
    const bool point = qu.lo == qu.hi && qv.lo == qv.hi;
    const double q[2] = {num::mid(qu), num::mid(qv)};

    for (std::size_t i = 0; i < cyc.links.size();) {
        const Link& L = cyc.links[i];
        if (num::disjoint(L.bu, qu) || num::disjoint(L.bv, qv)) { ++i; continue; }
        if (!L.arc) {
            // With exact ends the chord is the piece, and only a point on it is
            // undecided. With enclosed ends (or a query box) the piece lies
            // somewhere in the hull of its end boxes: the query must be
            // provably clear of that hull.
            if (point && L.exact_ends) {
                if (num::orient2d(L.p0, L.p1, q) == num::Sign::Zero) { res.on_boundary = true; return res; }
            } else if (!clear_of_straight(L, qu, qv)) {
                return res;
            }
            ++i;
            continue;
        }
        if (!L.arc->refinable || L.depth >= max_depth) return res;
        Link a, b;
        split_link(L, a, b);
        cyc.links[i] = a;
        cyc.links.insert(cyc.links.begin() + static_cast<std::ptrdiff_t>(i) + 1, b);
    }

    int wn = 0;
    for (const Link& L : cyc.links) {
        if (L.p0[1] <= q[1]) {
            if (L.p1[1] > q[1]) {
                const num::Sign s = num::orient2d(L.p0, L.p1, q);
                if (s == num::Sign::Zero) return res;
                if (s == num::Sign::Positive) ++wn;
            }
        } else if (L.p1[1] <= q[1]) {
            const num::Sign s = num::orient2d(L.p0, L.p1, q);
            if (s == num::Sign::Zero) return res;
            if (s == num::Sign::Negative) --wn;
        }
    }
    res.decided = true;
    res.wn = wn;
    return res;
}

/// Orientation of a loop that is a simple closed curve: +1 counter-clockwise,
/// -1 clockwise, 0 if it could not be decided. Found by locating any point
/// with a nonzero certified winding number.
inline int cycle_orientation(CycleGeom& cyc, int max_depth) {
    struct Chord { double mx, my, nx, ny, len; };
    std::vector<Chord> chords;
    double cx = 0.0, cy = 0.0;
    for (const Link& L : cyc.links) {
        const double dx = L.p1[0] - L.p0[0], dy = L.p1[1] - L.p0[1];
        const double len = std::sqrt(dx * dx + dy * dy);
        cx += L.p0[0];
        cy += L.p0[1];
        if (len > 0.0)
            chords.push_back({0.5 * (L.p0[0] + L.p1[0]), 0.5 * (L.p0[1] + L.p1[1]), -dy / len, dx / len, len});
    }
    if (cyc.links.empty()) return 0;
    cx /= static_cast<double>(cyc.links.size());
    cy /= static_cast<double>(cyc.links.size());

    auto probe = [&](double x, double y) {
        const WindingResult w = winding(cyc, Interval(x), Interval(y), max_depth);
        return w.decided ? w.wn : 0;
    };
    if (const int w = probe(cx, cy)) return w > 0 ? 1 : -1;
    const std::size_t stride = std::max<std::size_t>(1, chords.size() / 6);
    for (std::size_t k = 0; k < chords.size(); k += stride)
        for (double frac : {0.25, 0.05, 0.01, 0.002})
            for (double sgn : {1.0, -1.0}) {
                const Chord& c = chords[k];
                const double d = sgn * frac * c.len;
                if (const int w = probe(c.mx + d * c.nx, c.my + d * c.ny)) return w > 0 ? 1 : -1;
            }
    return 0;
}

} // namespace ak::op
