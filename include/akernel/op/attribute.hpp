// akernel — attributes: named data attached to a body, a face, an edge or a
// vertex (Parasolid's attribute system, docs/FEATURES.md).
//
// An attribute is plain data stored by value in the entity it is attached
// to, so copying a Brep copies its attributes and every operation that
// copies an entity keeps them. Each field has a declared behaviour under a
// rigid map, so that placing a body keeps its attributes meaningful:
//
//   ints, reals, text   unchanged
//   points              positions: mapped like the body (x' = M x + t)
//   directions          free vectors: mapped by the linear part only
//
// Under the kernel's exact maps (op/brep_transform.hpp) both are exact.
#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include <akernel/geom/vec.hpp>

namespace ak::op {

struct Attribute {
    std::string name; ///< what the attribute is ("colour", "material", "feature id", ...)
    std::vector<std::int64_t> ints;
    std::vector<double> reals;
    std::string text;
    std::vector<geom::Vec3d> points;     ///< positions, moved with the body
    std::vector<geom::Vec3d> directions; ///< free vectors, turned with the body
};

using Attributes = std::vector<Attribute>;

inline bool same_vec(const geom::Vec3d& a, const geom::Vec3d& b) { return a.x == b.x && a.y == b.y && a.z == b.z; }

/// Field by field, exact.
inline bool same_attribute(const Attribute& a, const Attribute& b) {
    if (a.name != b.name || a.ints != b.ints || a.reals != b.reals || a.text != b.text) return false;
    if (a.points.size() != b.points.size() || a.directions.size() != b.directions.size()) return false;
    for (std::size_t i = 0; i < a.points.size(); ++i)
        if (!same_vec(a.points[i], b.points[i])) return false;
    for (std::size_t i = 0; i < a.directions.size(); ++i)
        if (!same_vec(a.directions[i], b.directions[i])) return false;
    return true;
}

/// The attribute with this name, or null.
inline const Attribute* find_attribute(const Attributes& as, const std::string& name) {
    for (const Attribute& a : as)
        if (a.name == name) return &a;
    return nullptr;
}

/// Attach a, replacing an attribute of the same name: an entity holds at
/// most one attribute per name.
inline void set_attribute(Attributes& as, Attribute a) {
    for (Attribute& b : as)
        if (b.name == a.name) { b = std::move(a); return; }
    as.push_back(std::move(a));
}

/// Remove the attribute with this name; false if there was none.
inline bool remove_attribute(Attributes& as, const std::string& name) {
    for (std::size_t i = 0; i < as.size(); ++i)
        if (as[i].name == name) {
            as.erase(as.begin() + static_cast<std::ptrdiff_t>(i));
            return true;
        }
    return false;
}

/// Apply a rigid map to every attribute: points through `point`, directions
/// through `linear`.
template <class PointMap, class LinearMap>
void map_attributes(Attributes& as, const PointMap& point, const LinearMap& linear) {
    for (Attribute& a : as) {
        for (geom::Vec3d& p : a.points) p = point(p);
        for (geom::Vec3d& d : a.directions) d = linear(d);
    }
}

} // namespace ak::op
