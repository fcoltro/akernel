// akernel — small fixed-size vectors, generic over the scalar type so the
// same geometry code runs in doubles and in intervals.
#pragma once

#include <cmath>
#include <type_traits>

namespace ak::geom {

template <class T>
struct Vec3 {
    T x{}, y{}, z{};
};

template <class T>
struct Vec4 {
    T x{}, y{}, z{}, w{};
};

// The scalar operand is a non-deduced context, so Vec3<Interval> * 2.0 works.
template <class T> using Scalar = std::type_identity_t<T>;

template <class T> Vec3<T> operator+(const Vec3<T>& a, const Vec3<T>& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
template <class T> Vec3<T> operator-(const Vec3<T>& a, const Vec3<T>& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
template <class T> Vec3<T> operator-(const Vec3<T>& a) { return {-a.x, -a.y, -a.z}; }
template <class T> Vec3<T> operator*(const Vec3<T>& a, const Scalar<T>& s) { return {a.x * s, a.y * s, a.z * s}; }
template <class T> Vec3<T> operator*(const Scalar<T>& s, const Vec3<T>& a) { return {a.x * s, a.y * s, a.z * s}; }
template <class T> Vec3<T> operator/(const Vec3<T>& a, const Scalar<T>& s) { return {a.x / s, a.y / s, a.z / s}; }
template <class T> Vec3<T>& operator+=(Vec3<T>& a, const Vec3<T>& b) { a = a + b; return a; }
template <class T> Vec3<T>& operator-=(Vec3<T>& a, const Vec3<T>& b) { a = a - b; return a; }
template <class T> T dot(const Vec3<T>& a, const Vec3<T>& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
template <class T> Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

template <class T> Vec4<T> operator+(const Vec4<T>& a, const Vec4<T>& b) { return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w}; }
template <class T> Vec4<T> operator-(const Vec4<T>& a, const Vec4<T>& b) { return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w}; }
template <class T> Vec4<T> operator*(const Vec4<T>& a, const Scalar<T>& s) { return {a.x * s, a.y * s, a.z * s, a.w * s}; }
template <class T> Vec4<T> operator*(const Scalar<T>& s, const Vec4<T>& a) { return {a.x * s, a.y * s, a.z * s, a.w * s}; }

template <class T> Vec3<T> xyz(const Vec4<T>& a) { return {a.x, a.y, a.z}; }

using Vec3d = Vec3<double>;
using Vec4d = Vec4<double>;

inline double norm(const Vec3d& a) { return std::sqrt(dot(a, a)); }
inline Vec3d normalized(const Vec3d& a) { return a / norm(a); }

} // namespace ak::geom
