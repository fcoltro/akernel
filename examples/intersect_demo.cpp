// Intersect two freeform patches and write the result as a Wavefront OBJ:
// both surfaces as quad grids and the certified intersection curve as a
// polyline. Open the file in any mesh viewer.
//
//   ak_demo [output.obj]
#include <akernel/geom/primitives.hpp>
#include <akernel/solve/ssi.hpp>

#include <cstdio>
#include <fstream>
#include <iomanip>

using namespace ak;
using geom::Vec3d;

int main(int argc, char** argv) {
    const char* path = argc > 1 ? argv[1] : "akernel_demo.obj";

    const geom::PolyPatchD a = geom::fit_patch(
        [](double u, double v) { return Vec3d{u, v, 0.35 * u * u * u - 0.5 * u * v + 0.3 * v * v + 0.1}; }, 3, 3);
    const geom::PolyPatchD b = geom::fit_patch(
        [](double s, double t) {
            return Vec3d{0.05 + 0.9 * s + 0.1 * t * t, 0.95 * t + 0.02, 0.7 * s - 0.25 * t * s * s + 0.05 * t - 0.05};
        }, 3, 3);

    const solve::SsiSystem sys(a, b);
    const solve::SsiResult res = solve::intersect(sys);

    std::ofstream f(path);
    if (!f) { std::printf("cannot write %s\n", path); return 1; }
    f << std::setprecision(17) << "# akernel intersection demo\n";

    int next = 1; // OBJ vertex indices start at 1
    const int grid = 40;
    const geom::BezierSurface surf[2] = {geom::BezierSurface(a), geom::BezierSurface(b)};
    for (int s = 0; s < 2; ++s) {
        f << "o surface_" << (s == 0 ? 'a' : 'b') << "\n";
        for (int i = 0; i <= grid; ++i)
            for (int j = 0; j <= grid; ++j) {
                const Vec3d p = surf[s].point(static_cast<double>(i) / grid, static_cast<double>(j) / grid);
                f << "v " << p.x << ' ' << p.y << ' ' << p.z << "\n";
            }
        for (int i = 0; i < grid; ++i)
            for (int j = 0; j < grid; ++j) {
                const int v00 = next + i * (grid + 1) + j;
                f << "f " << v00 << ' ' << v00 + grid + 1 << ' ' << v00 + grid + 2 << ' ' << v00 + 1 << "\n";
            }
        next += (grid + 1) * (grid + 1);
    }

    std::size_t curve_id = 0, points = 0;
    double worst = 0.0;
    for (const auto& c : res.curves) {
        f << "o intersection_" << curve_id++ << "\n";
        const int first = next;
        for (const auto& seg : c.segs)
            for (double t : {0.0, 0.25, 0.5, 0.75}) {
                solve::Box4 enc;
                if (!solve::curve_point(sys, seg, t, enc)) continue;
                const auto pa = sys.point_a(enc), pb = sys.point_b(enc);
                worst = std::max(worst, std::sqrt((sqr(pa.x - pb.x) + sqr(pa.y - pb.y) + sqr(pa.z - pb.z)).hi));
                f << "v " << mid(pa.x) << ' ' << mid(pa.y) << ' ' << mid(pa.z) << "\n";
                ++next;
                ++points;
            }
        f << "l";
        for (int v = first; v < next; ++v) f << ' ' << v;
        if (c.closed) f << ' ' << first;
        f << "\n";
    }
    f.close();

    std::printf("wrote %s\n", path);
    std::printf("curves: %zu, curve points: %zu, proven complete: %s\n", res.curves.size(), points,
                res.complete ? "yes" : "no");
    std::printf("certified bound on the distance between the two surfaces at every curve point: %.1e\n", worst);
    return 0;
}
