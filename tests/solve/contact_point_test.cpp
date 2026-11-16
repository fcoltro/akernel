// Tangent contacts classified with certificates (solve/contact_point.hpp):
// a paraboloid cap against planes at, just above and just below its apex,
// a saddle against its tangent plane, and the same in rational form.
#include <akernel/solve/contact_point.hpp>

#include <cmath>
#include <cstdio>

#include "check.hpp"

using namespace ak;
using geom::PolyPatchD;
using solve::ContactKind;

namespace {

/// z = 1 - (x^2 + y^2) / 2 over x, y in [-1, 1] (u, v in [0, 1]); with
/// `saddle`, z = (x^2 - y^2) / 2. Apex at u = v = 1/2.
PolyPatchD quadric_cap(bool saddle, double weight = 1.0) {
    PolyPatchD P(2, 2);
    const double a[3] = {1.0, -1.0, 1.0}; // Bernstein form of (2u - 1)^2
    for (int i = 0; i <= 2; ++i)
        for (int j = 0; j <= 2; ++j) {
            const double z = saddle ? 0.5 * (a[i] - a[j]) : 1.0 - 0.5 * (a[i] + a[j]);
            P.at(i, j) = geom::Vec4d{-1.0 + i, -1.0 + j, z, 1.0} * weight;
        }
    return P;
}

/// The plane z = c over x, y in [-0.75, 1.25]: the apex is at s = t = 3/8.
PolyPatchD plane(double c) {
    PolyPatchD P(1, 1);
    for (int i = 0; i <= 1; ++i)
        for (int j = 0; j <= 1; ++j) P.at(i, j) = {-0.75 + 2.0 * i, -0.75 + 2.0 * j, c, 1.0};
    return P;
}

const std::array<double, 5> kSeed = {0.47, 0.52, 0.40, 0.36, 0.01};

void report(const char* name, const solve::ContactResult& r) {
    std::printf("  %s: %s", name, solve::to_string(r.kind));
    if (r.kind != ContactKind::NotFound)
        std::printf(", mu in [%.3e, %.3e], relative curvature %s", r.enc[4].lo, r.enc[4].hi,
                    r.d_sign == 1 ? "positive definite" : r.d_sign == -1 ? "negative definite" : r.d_sign == 2 ? "indefinite" : "undecided");
    if (r.exact) std::printf(", exactly at (u, v, s, t) = (%g, %g, %g, %g)", r.point[0], r.point[1], r.point[2], r.point[3]);
    if (!r.why.empty()) std::printf(" (%s)", r.why.c_str());
    std::printf("\n");
}

} // namespace

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);

    // 1. The plane through the apex touches the cap there, and nowhere else
    //    nearby: an isolated tangent point, proven exactly.
    {
        const auto r = solve::classify_contact(quadric_cap(false), plane(1.0), kSeed);
        report("cap against the plane through its apex", r);
        AK_CHECK(r.kind == ContactKind::Touch);
        AK_CHECK(r.exact && r.point[0] == 0.5 && r.point[1] == 0.5 && r.point[2] == 0.375 && r.point[3] == 0.375);
        AK_CHECK(r.d_sign == 1);
    }
    // 2. Raised by 2^-20: the facing points are apart, the plane stays above.
    {
        const auto r = solve::classify_contact(quadric_cap(false), plane(1.0 + std::ldexp(1.0, -20)), kSeed);
        report("plane 2^-20 above the apex", r);
        AK_CHECK(r.kind == ContactKind::Separated);
        AK_CHECK(r.enc[4].lo > 0.0);
    }
    // 3. Lowered by 2^-20: the plane cuts a small circle out of the cap.
    {
        const auto r = solve::classify_contact(quadric_cap(false), plane(1.0 - std::ldexp(1.0, -20)), kSeed);
        report("plane 2^-20 below the apex", r);
        AK_CHECK(r.kind == ContactKind::Loop);
        AK_CHECK(r.enc[4].hi < 0.0);
    }
    // 4. Raised by 2^-44 (about 6e-14): still decided.
    {
        const auto r = solve::classify_contact(quadric_cap(false), plane(1.0 + std::ldexp(1.0, -44)), kSeed);
        report("plane 2^-44 above the apex", r);
        AK_CHECK(r.kind == ContactKind::Separated);
    }
    // 5. A saddle against its tangent plane: they touch exactly at the
    //    saddle point, where two intersection curves cross.
    {
        const auto r = solve::classify_contact(quadric_cap(true), plane(0.0), kSeed);
        report("saddle against its tangent plane", r);
        AK_CHECK(r.kind == ContactKind::SingularCrossing);
        AK_CHECK(r.exact && r.d_sign == 2);
    }
    // 6. The cap with every weight 2 (the same surface, as a rational patch):
    //    the same answers through the homogeneous computation.
    {
        const auto t = solve::classify_contact(quadric_cap(false, 2.0), plane(1.0), kSeed);
        report("rational form of the cap, touching", t);
        AK_CHECK(t.kind == ContactKind::Touch && t.exact);
        const auto l = solve::classify_contact(quadric_cap(false, 2.0), plane(1.0 - std::ldexp(1.0, -20)), kSeed);
        report("rational form of the cap, plane below", l);
        AK_CHECK(l.kind == ContactKind::Loop);
    }
    // 6b. A genuinely rational cap: weights (1, 2, 1) x (1, 2, 1), symmetric,
    //     so the apex stays at u = v = 1/2, x = y = 0. With z coefficients
    //     1 - 0.375 (a_i + a_j) its height there is 1 + 0.75 / 3 = 5/4
    //     (the weighted mean of a is -1/3), so the plane z = 5/4 touches it.
    {
        PolyPatchD P(2, 2);
        const double a[3] = {1.0, -1.0, 1.0}, w[3] = {1.0, 2.0, 1.0};
        for (int i = 0; i <= 2; ++i)
            for (int j = 0; j <= 2; ++j)
                P.at(i, j) = geom::Vec4d{-1.0 + i, -1.0 + j, 1.0 - 0.375 * (a[i] + a[j]), 1.0} * (w[i] * w[j]);
        const auto t = solve::classify_contact(P, plane(1.25), kSeed);
        report("rational cap (weights 1, 2, 1), plane through its apex", t);
        AK_CHECK(t.kind == ContactKind::Touch && t.exact && t.point[0] == 0.5 && t.point[1] == 0.5);
        const auto s = solve::classify_contact(P, plane(1.25 + std::ldexp(1.0, -30)), kSeed);
        report("rational cap, plane 2^-30 above", s);
        AK_CHECK(s.kind == ContactKind::Separated);
    }
    // 7. A cap whose apex is not at dyadic parameters (shifted by 1/3 in x):
    //    the contact is still proven unique in a box, but its exactness
    //    cannot be shown with dyadic candidates.
    {
        PolyPatchD P = quadric_cap(false);
        for (auto& c : P.c) c.x += 1.0 / 3.0;
        const auto r = solve::classify_contact(P, plane(1.0), {0.5, 0.5, 0.55, 0.37, 0.0});
        report("cap with its apex at x = 1/3", r);
        AK_CHECK(r.kind == ContactKind::Undecided);
    }
    return akt::finish("contact_point_test");
}
