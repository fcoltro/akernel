// akernel — a static bounding-volume hierarchy over certified boxes.
//
// The broad phase of the Boolean asks the same question many times: which of
// these n boxes may meet that one? Each box encloses an entity (a face, an
// edge, a vertex) for certain, so two entities whose boxes are disjoint
// cannot meet, and skipping them changes no answer. The index answers the
// question in about log n + (number of hits) box tests instead of n.
//
// Built once, top-down: the boxes of a node are split at the median of their
// centres along the axis where the centres spread most. Queries report ids
// in increasing order, so code that relied on scanning ids in order keeps
// the same order (and therefore the same results).
#pragma once

#include <algorithm>
#include <vector>

#include <akernel/geom/vec.hpp>
#include <akernel/num/interval.hpp>

namespace ak::geom {

class BoxIndex {
public:
    using Box = Vec3<num::Interval>;

    BoxIndex() = default;
    explicit BoxIndex(std::vector<Box> boxes) : boxes_(std::move(boxes)) {
        const int n = static_cast<int>(boxes_.size());
        order_.resize(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) order_[static_cast<std::size_t>(i)] = i;
        if (n > 0) build(0, n);
    }

    std::size_t size() const { return boxes_.size(); }
    const Box& box(int id) const { return boxes_[static_cast<std::size_t>(id)]; }

    /// Same rule as the rest of the kernel: closed boxes, disjoint only when
    /// separated along some axis.
    static bool disjoint(const Box& a, const Box& b) {
        return num::disjoint(a.x, b.x) || num::disjoint(a.y, b.y) || num::disjoint(a.z, b.z);
    }

    /// Ids of all boxes not disjoint from q, in increasing order.
    std::vector<int> query(const Box& q) const {
        std::vector<int> hits;
        if (nodes_.empty()) return hits;
        std::vector<int> stack{0};
        while (!stack.empty()) {
            const Node& nd = nodes_[static_cast<std::size_t>(stack.back())];
            stack.pop_back();
            if (disjoint(nd.box, q)) continue;
            if (nd.left < 0) {
                for (int i = nd.first; i < nd.first + nd.count; ++i) {
                    const int id = order_[static_cast<std::size_t>(i)];
                    if (!disjoint(boxes_[static_cast<std::size_t>(id)], q)) hits.push_back(id);
                }
                continue;
            }
            stack.push_back(nd.left);
            stack.push_back(nd.right);
        }
        std::sort(hits.begin(), hits.end());
        return hits;
    }

private:
    struct Node {
        Box box;
        int left = -1, right = -1; ///< children, or -1 for a leaf
        int first = 0, count = 0;  ///< a leaf's range in order_
    };
    static constexpr int kLeaf = 4;

    std::vector<Box> boxes_;
    std::vector<int> order_;
    std::vector<Node> nodes_;

    static num::Interval hull(const num::Interval& a, const num::Interval& b) {
        return {std::min(a.lo, b.lo), std::max(a.hi, b.hi)};
    }
    static double centre(const Box& b, int axis) {
        const num::Interval& c = axis == 0 ? b.x : axis == 1 ? b.y : b.z;
        return 0.5 * c.lo + 0.5 * c.hi;
    }

    int build(int first, int count) {
        const int id = static_cast<int>(nodes_.size());
        nodes_.push_back({});
        Box b = boxes_[static_cast<std::size_t>(order_[static_cast<std::size_t>(first)])];
        double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
        for (int i = first; i < first + count; ++i) {
            const Box& c = boxes_[static_cast<std::size_t>(order_[static_cast<std::size_t>(i)])];
            b = {hull(b.x, c.x), hull(b.y, c.y), hull(b.z, c.z)};
            for (int a = 0; a < 3; ++a) {
                lo[a] = std::min(lo[a], centre(c, a));
                hi[a] = std::max(hi[a], centre(c, a));
            }
        }
        nodes_[static_cast<std::size_t>(id)].box = b;
        if (count <= kLeaf) {
            nodes_[static_cast<std::size_t>(id)].first = first;
            nodes_[static_cast<std::size_t>(id)].count = count;
            return id;
        }
        int axis = 0;
        for (int a = 1; a < 3; ++a)
            if (hi[a] - lo[a] > hi[axis] - lo[axis]) axis = a;
        const int half = count / 2;
        auto begin = order_.begin() + first;
        std::nth_element(begin, begin + half, begin + count, [&](int p, int q) {
            return centre(boxes_[static_cast<std::size_t>(p)], axis) < centre(boxes_[static_cast<std::size_t>(q)], axis);
        });
        const int l = build(first, half);
        const int r = build(first + half, count - half);
        nodes_[static_cast<std::size_t>(id)].left = l;
        nodes_[static_cast<std::size_t>(id)].right = r;
        return id;
    }
};

} // namespace ak::geom
