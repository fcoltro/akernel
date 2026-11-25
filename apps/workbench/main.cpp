// ak_workbench — a minimal desktop workbench for akernel.
//
// A 3D viewport (orbit, pan, zoom, click to select) and panels of buttons
// that call the kernel's own operations: primitives, extrude and revolve
// (with blended edges), sweep and loft surfaces, the freeform fillet,
// Booleans, plane cuts and exact placement. Each result shows what the
// kernel certified or measured, or the kernel's reason for refusing.
//
// Kernel calls run on a worker thread so the window stays responsive while
// a Boolean takes seconds. The kernel calls themselves are in kernel_ops.hpp.
//
// Mouse: left drag orbits, right or middle drag pans, wheel zooms. Click a
// body to make it A, Shift+click to make it B. F fits the view.
//
// Command line (for scripted checks): --demo builds a small scene at start;
// --shot file.ppm renders it once, saves the window as a PPM image and exits.
#include "kernel_ops.hpp"

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl2.h>

#include <GLFW/glfw3.h>

#include <atomic>
#include <cmath>
#include <deque>
#include <exception>
#include <fstream>
#include <functional>
#include <mutex>
#include <numbers>
#include <thread>

#ifndef GL_MULTISAMPLE
#define GL_MULTISAMPLE 0x809D
#endif

namespace {

using wb::Body;
using wb::Outcome;
using wb::Vec3d;

// ------------------------------------------------------------------ worker

/// Runs one kernel call at a time off the UI thread.
class Worker {
public:
    ~Worker() {
        if (th_.joinable()) th_.join();
    }
    bool busy() const { return busy_; }
    const std::string& label() const { return label_; }
    double elapsed() const {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0_).count();
    }
    bool start(std::string label, std::function<Outcome()> f) {
        if (busy_) return false;
        if (th_.joinable()) th_.join();
        busy_ = true;
        label_ = std::move(label);
        t0_ = std::chrono::steady_clock::now();
        th_ = std::thread([this, job = std::move(f)] {
            Outcome o;
            try {
                o = wb::timed(job);
            } catch (const std::exception& e) {
                o = wb::refused(std::string("exception: ") + e.what());
            }
            {
                std::lock_guard<std::mutex> lk(m_);
                done_ = std::move(o);
            }
            busy_ = false;
        });
        return true;
    }
    std::optional<Outcome> poll() {
        std::lock_guard<std::mutex> lk(m_);
        std::optional<Outcome> r = std::move(done_);
        done_.reset();
        return r;
    }

private:
    std::thread th_;
    std::atomic<bool> busy_{false};
    std::mutex m_;
    std::optional<Outcome> done_;
    std::string label_;
    std::chrono::steady_clock::time_point t0_;
};

// ------------------------------------------------------------------ camera

struct Camera {
    Vec3d target{0.5, 0.5, 0.25};
    double yaw = -0.9, pitch = 0.5, dist = 4.0;
    double fovy = 0.6; ///< radians

    Vec3d eye() const {
        return target + Vec3d{std::cos(pitch) * std::cos(yaw), std::cos(pitch) * std::sin(yaw), std::sin(pitch)} * dist;
    }
    /// Orthonormal frame: forward (eye to target), right, up.
    void frame(Vec3d& f, Vec3d& s, Vec3d& u) const {
        f = target - eye();
        f = f / norm(f);
        s = cross(f, Vec3d{0, 0, 1});
        s = s / norm(s);
        u = cross(s, f);
    }
    void fit(const std::array<float, 3>& lo, const std::array<float, 3>& hi) {
        target = {0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])};
        const double r = 0.5 * norm(Vec3d{double(hi[0]) - lo[0], double(hi[1]) - lo[1], double(hi[2]) - lo[2]});
        dist = std::max(r, 1e-3) / std::sin(0.5 * fovy) * 1.1;
    }
    void apply(int w, int h) const {
        const double aspect = h > 0 ? double(w) / h : 1.0;
        const double zn = dist * 0.01, zf = dist * 100.0;
        const double th = std::tan(0.5 * fovy) * zn;
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glFrustum(-th * aspect, th * aspect, -th, th, zn, zf);
        glMatrixMode(GL_MODELVIEW);
        Vec3d f, s, u;
        frame(f, s, u);
        const Vec3d e = eye();
        const double m[16] = {s.x, u.x, -f.x, 0, s.y, u.y, -f.y, 0, s.z, u.z, -f.z, 0,
                              -dot(s, e), -dot(u, e), dot(f, e), 1};
        glLoadMatrixd(m);
    }
    /// The ray through pixel (x, y) of a w x h viewport (y down).
    void ray(double x, double y, int w, int h, Vec3d& o, Vec3d& d) const {
        Vec3d f, s, u;
        frame(f, s, u);
        const double t = std::tan(0.5 * fovy), aspect = h > 0 ? double(w) / h : 1.0;
        const double px = (2.0 * x / w - 1.0) * t * aspect, py = (1.0 - 2.0 * y / h) * t;
        o = eye();
        d = f + s * px + u * py;
        d = d / norm(d);
    }
};

// ------------------------------------------------------------------- app

struct LogLine {
    bool ok;
    std::string text;
};

struct App {
    std::vector<Body> bodies;
    int selA = -1, selB = -1;
    std::deque<LogLine> log;
    Worker worker;
    Camera cam;
    bool hide_operands = true;
    bool show_edges = true, show_grid = true;
    std::vector<int> pending_operands; ///< bodies to hide when the running job succeeds

    // Tool parameters (dyadic defaults, so the kernel's controls stay exact).
    Vec3d origin{0, 0, 0}, size{1, 1, 0.5};
    double bump = 0.25, round_r = 0.125, amp = 0.0625;
    int spans_m = 2, spans_n = 2;
    double ex_w = 1, ex_h = 0.5, ex_l = 1, ex_r = 0.125, ex_d2 = 0.125, step_a1 = 0.5, step_b1 = 0.25;
    int ex_corner = 1;
    double rv_r0 = 0.25, rv_r1 = 0.5, rv_z0 = 0, rv_z1 = 0.5, rv_r = 0.125, rv_d2 = 0.125;
    int rv_quarters = 4, rv_edge = 0, rv_corner = 0;
    double sw_r = 0.25;
    Vec3d sw_c1{0, 0.5, 0.5}, sw_c2{0.5, 0.5, 1}, sw_end{1, 0, 1.5};
    double lf_r0 = 0.5, lf_r1 = 0.25, lf_r2 = 0.375, lf_h = 1;
    double ff_a = 1, ff_L = 2, ff_h = 0.5, ff_c = 0.25, ff_r = 0.125, ff_cut = 1;
    Vec3d cut_n{0, 0, 1};
    double cut_d = 0.25;
    Vec3d move_t{0.5, 0, 0};

    void say(bool ok, std::string s) {
        log.push_front({ok, std::move(s)});
        if (log.size() > 200) log.pop_back();
    }

    void run(std::string label, std::function<Outcome()> f, std::vector<int> operands = {}) {
        if (worker.start(label, std::move(f))) pending_operands = std::move(operands);
        else say(false, "Busy: wait for \"" + worker.label() + "\" to finish");
    }

    void add(Body b) {
        static const std::array<float, 3> palette[] = {{0.62f, 0.70f, 0.82f}, {0.78f, 0.72f, 0.60f},
                                                       {0.64f, 0.78f, 0.66f}, {0.80f, 0.64f, 0.70f},
                                                       {0.70f, 0.66f, 0.82f}, {0.62f, 0.78f, 0.80f}};
        if (b.solid) b.color = palette[bodies.size() % 6];
        bodies.push_back(std::move(b));
        selA = static_cast<int>(bodies.size()) - 1;
        if (selB == selA) selB = -1;
    }

    void poll() {
        std::optional<Outcome> o = worker.poll();
        if (!o) return;
        if (o->ok && o->body) {
            for (int i : pending_operands)
                if (hide_operands && i >= 0 && i < static_cast<int>(bodies.size())) bodies[static_cast<std::size_t>(i)].visible = false;
            say(true, wb::fmt("%s  (%.2f s)", o->message.c_str(), o->seconds));
            const bool first = bodies.empty();
            add(std::move(*o->body));
            if (first) fit_all();
        } else {
            say(false, wb::fmt("Refused: %s  (%.2f s)", o->message.c_str(), o->seconds));
        }
        pending_operands.clear();
    }

    void fit_all() {
        std::array<float, 3> lo{1e30f, 1e30f, 1e30f}, hi{-1e30f, -1e30f, -1e30f};
        bool any = false;
        for (const Body& b : bodies)
            if (b.visible && !b.mesh.empty()) {
                any = true;
                for (std::size_t k = 0; k < 3; ++k) {
                    lo[k] = std::min(lo[k], b.mesh.lo[k]);
                    hi[k] = std::max(hi[k], b.mesh.hi[k]);
                }
            }
        if (any) cam.fit(lo, hi);
    }

    const Body* body(int i) const {
        return i >= 0 && i < static_cast<int>(bodies.size()) ? &bodies[static_cast<std::size_t>(i)] : nullptr;
    }

    /// Nearest visible body hit by the ray, or -1.
    int pick(const Vec3d& o, const Vec3d& d) const {
        double best = 1e300;
        int hit = -1;
        for (std::size_t bi = 0; bi < bodies.size(); ++bi) {
            const Body& b = bodies[bi];
            if (!b.visible) continue;
            const auto& t = b.mesh.tri;
            for (std::size_t k = 0; k + 18 <= t.size(); k += 18) {
                const Vec3d p0{t[k], t[k + 1], t[k + 2]}, p1{t[k + 6], t[k + 7], t[k + 8]}, p2{t[k + 12], t[k + 13], t[k + 14]};
                const Vec3d e1 = p1 - p0, e2 = p2 - p0, pv = cross(d, e2);
                const double det = dot(e1, pv);
                if (std::fabs(det) < 1e-14) continue;
                const Vec3d tv = o - p0;
                const double u = dot(tv, pv) / det;
                if (u < 0 || u > 1) continue;
                const Vec3d qv = cross(tv, e1);
                const double v = dot(d, qv) / det;
                if (v < 0 || u + v > 1) continue;
                const double s = dot(e2, qv) / det;
                if (s > 0 && s < best) {
                    best = s;
                    hit = static_cast<int>(bi);
                }
            }
        }
        return hit;
    }

    void demo() {
        auto take = [&](Outcome o) {
            if (o.ok && o.body) add(std::move(*o.body));
            else say(false, "Demo: " + o.message);
        };
        // A plate with a ring-shaped groove cut by a Boolean, a filleted
        // extrusion, a rounded box, the freeform fillet and a loft.
        Outcome plate = wb::box({0, 0, 0}, {1, 1, 0.5});
        Outcome ring = wb::revolve_rect(0.25, 0.375, 0.25, 0.75, 4, 0, wb::Corner::None, 0, 0);
        if (plate.ok && ring.ok) {
            Outcome moved = wb::place(*ring.body, ak::op::SignedPermutation::translation({0.5, 0.5, 0}), "");
            if (moved.ok) {
                Outcome diff = wb::boolean(*plate.body, *moved.body, 1);
                take(diff.ok ? std::move(diff) : std::move(plate));
            }
        }
        take(wb::extrude_rect({1.5, 0, 0}, 0.5, 0.5, 1, wb::Corner::Fillet, 0.125, 0));
        take(wb::rounded_box({2.5, 0, 0}, {0.75, 0.75, 0.75}, 0.125));
        {
            Outcome rv = wb::revolve_rect(0.25, 0.5, 0, 0.375, 3, 0, wb::Corner::Fillet, 0.125, 0);
            if (rv.ok) take(wb::place(*rv.body, ak::op::SignedPermutation::translation({1, 2, 0}), ""));
            else take(std::move(rv));
        }
        take(wb::loft_circles({3.25, 1.25, 0}, 0.25, 0.125, 0.25, 0.75));
        {
            Outcome ff = wb::freeform_fillet(1, 2, 0.5, 0.25, 0.125, 1);
            if (ff.ok) take(wb::place(*ff.body, ak::op::SignedPermutation::translation({-1.5, 0, 0}), ""));
            else take(std::move(ff));
        }
        fit_all();
    }
};

// --------------------------------------------------------------- drawing

void draw_grid() {
    glDisable(GL_LIGHTING);
    glLineWidth(1.0f);
    glBegin(GL_LINES);
    glColor3f(0.32f, 0.33f, 0.36f);
    for (int i = -10; i <= 10; ++i) {
        const float t = 0.25f * static_cast<float>(i);
        glVertex3f(t, -2.5f, 0);
        glVertex3f(t, 2.5f, 0);
        glVertex3f(-2.5f, t, 0);
        glVertex3f(2.5f, t, 0);
    }
    glEnd();
    glLineWidth(2.0f);
    glBegin(GL_LINES);
    glColor3f(0.85f, 0.3f, 0.3f);
    glVertex3f(0, 0, 0);
    glVertex3f(0.5f, 0, 0);
    glColor3f(0.3f, 0.8f, 0.3f);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0.5f, 0);
    glColor3f(0.35f, 0.5f, 0.95f);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, 0.5f);
    glEnd();
}

void draw_body(const Body& b, int role, bool edges) {
    if (b.mesh.empty()) return;
    const std::array<float, 3> tint = role == 1 ? std::array<float, 3>{0.35f, 0.6f, 1.0f}
                                    : role == 2 ? std::array<float, 3>{1.0f, 0.6f, 0.25f}
                                                : b.color;
    glEnable(GL_LIGHTING);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0f, 1.0f);
    glColor3f(tint[0], tint[1], tint[2]);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3, GL_FLOAT, 6 * sizeof(float), b.mesh.tri.data());
    glNormalPointer(GL_FLOAT, 6 * sizeof(float), b.mesh.tri.data() + 3);
    glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(b.mesh.tri.size() / 6));
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisable(GL_POLYGON_OFFSET_FILL);
    if (edges && !b.mesh.lines.empty()) {
        glDisable(GL_LIGHTING);
        glLineWidth(role ? 2.0f : 1.2f);
        if (role) glColor3f(tint[0] * 0.5f, tint[1] * 0.5f, tint[2] * 0.5f);
        else glColor3f(0.08f, 0.08f, 0.1f);
        glVertexPointer(3, GL_FLOAT, 0, b.mesh.lines.data());
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(b.mesh.lines.size() / 3));
    }
    glDisableClientState(GL_VERTEX_ARRAY);
}

void draw_scene(const App& app, int w, int h) {
    glViewport(0, 0, w, h);
    glClearColor(0.16f, 0.17f, 0.19f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
    // Head lights: set with the identity model-view so they follow the camera.
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    const GLfloat l0[4] = {0.3f, 0.6f, 1.0f, 0.0f}, l1[4] = {-0.5f, -0.4f, 0.3f, 0.0f};
    const GLfloat d0[4] = {0.8f, 0.8f, 0.8f, 1}, d1[4] = {0.3f, 0.3f, 0.32f, 1}, amb[4] = {0.22f, 0.22f, 0.24f, 1};
    const GLfloat spec[4] = {0.25f, 0.25f, 0.25f, 1}, none[4] = {0, 0, 0, 1};
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glLightfv(GL_LIGHT0, GL_POSITION, l0);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, d0);
    glLightfv(GL_LIGHT0, GL_SPECULAR, spec);
    glLightfv(GL_LIGHT1, GL_POSITION, l1);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, d1);
    glLightfv(GL_LIGHT1, GL_SPECULAR, none);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, amb);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 32.0f);
    glEnable(GL_NORMALIZE);
    glShadeModel(GL_SMOOTH);
    app.cam.apply(w, h);
    if (app.show_grid) draw_grid();
    for (std::size_t i = 0; i < app.bodies.size(); ++i) {
        const Body& b = app.bodies[i];
        if (!b.visible) continue;
        const int idx = static_cast<int>(i);
        draw_body(b, idx == app.selA ? 1 : idx == app.selB ? 2 : 0, app.show_edges);
    }
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
}

// ------------------------------------------------------------------- UI

bool drag3(const char* label, Vec3d& v, double speed = 0.0625) {
    double a[3] = {v.x, v.y, v.z};
    const bool ch = ImGui::DragScalarN(label, ImGuiDataType_Double, a, 3, static_cast<float>(speed), nullptr, nullptr, "%.4g");
    v = {a[0], a[1], a[2]};
    return ch;
}

bool drag1(const char* label, double& v, double speed = 0.0625, double lo = -1e9, double hi = 1e9) {
    return ImGui::DragScalar(label, ImGuiDataType_Double, &v, static_cast<float>(speed), &lo, &hi, "%.4g");
}

void help(const char* text) {
    ImGui::SameLine();
    ImGui::TextDisabled("(?)");
    if (ImGui::BeginItemTooltip()) {
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 28.0f);
        ImGui::TextUnformatted(text);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}

const char* kCorners[] = {"None", "Fillet", "Chamfer"};

void corner_inputs(const char* id, int& corner, double& r, double& d2) {
    ImGui::PushID(id);
    ImGui::Combo("Corner", &corner, kCorners, 3);
    if (corner == 1) drag1("Radius", r, 0.015625, 0.0, 1e9);
    if (corner == 2) {
        drag1("Setback A", r, 0.015625, 0.0, 1e9);
        drag1("Setback B", d2, 0.015625, 0.0, 1e9);
    }
    ImGui::PopID();
}

void tools_panel(App& app, float height) {
    ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(330, height), ImGuiCond_Always);
    ImGui::Begin("Tools", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);
    ImGui::PushItemWidth(170);

    if (ImGui::CollapsingHeader("Primitives", ImGuiTreeNodeFlags_DefaultOpen)) {
        drag3("Origin", app.origin);
        drag3("Size", app.size);
        if (ImGui::Button("Box")) { const Vec3d o = app.origin, s = app.size; app.run("Box", [=] { return wb::box(o, s); }); }
        ImGui::SameLine();
        if (ImGui::Button("Bump box")) { const Vec3d o = app.origin, s = app.size; const double b = app.bump; app.run("Bump box", [=] { return wb::bump_box(o, s, b); }); }
        ImGui::SameLine();
        if (ImGui::Button("Rounded box")) { const Vec3d o = app.origin, s = app.size; const double r = app.round_r; app.run("Rounded box", [=] { return wb::rounded_box(o, s, r); }); }
        drag1("Bump height", app.bump, 0.015625);
        drag1("Round radius", app.round_r, 0.015625, 0.0, 1e9);
        if (ImGui::Button("Spline box")) { const Vec3d o = app.origin, s = app.size; const int m = app.spans_m, n = app.spans_n; const double a = app.amp; app.run("Spline box", [=] { return wb::spline_box(o, s, m, n, a); }); }
        help("A box whose top is a quadratic B-spline with m x n spans (one face per Bezier piece).");
        ImGui::SliderInt("Spans u", &app.spans_m, 1, 6);
        ImGui::SliderInt("Spans v", &app.spans_n, 1, 6);
        drag1("Wave height", app.amp, 0.015625);
    }

    if (ImGui::CollapsingHeader("Extrude")) {
        ImGui::TextWrapped("Rectangle W x H at Origin, extruded by L along z. The corner (W, H) can be rounded or chamfered.");
        drag1("W", app.ex_w, 0.0625, 0.0, 1e9);
        drag1("H", app.ex_h, 0.0625, 0.0, 1e9);
        drag1("L", app.ex_l, 0.0625, 0.0, 1e9);
        corner_inputs("ex", app.ex_corner, app.ex_r, app.ex_d2);
        if (ImGui::Button("Extrude")) {
            const Vec3d o = app.origin; const double W = app.ex_w, H = app.ex_h, L = app.ex_l, r = app.ex_r, d2 = app.ex_d2;
            const auto c = static_cast<wb::Corner>(app.ex_corner);
            app.run("Extrude", [=] { return wb::extrude_rect(o, W, H, L, c, r, d2); });
        }
        ImGui::Separator();
        drag1("Step a1", app.step_a1, 0.0625, 0.0, 1e9);
        drag1("Step b1", app.step_b1, 0.0625, 0.0, 1e9);
        if (ImGui::Button("Extrude step (inside edge)")) {
            const Vec3d o = app.origin; const double W = app.ex_w, H = app.ex_h, L = app.ex_l, a1 = app.step_a1, b1 = app.step_b1, r = app.ex_r, d2 = app.ex_d2;
            const auto c = static_cast<wb::Corner>(app.ex_corner);
            app.run("Step", [=] { return wb::extrude_step(o, W, H, L, a1, b1, c, r, d2); });
        }
        help("An L-shaped profile (W x H without the notch x > a1, y > b1). Its concave edge gets the corner fillet or chamfer: material is added.");
    }

    if (ImGui::CollapsingHeader("Revolve")) {
        ImGui::TextWrapped("Rectangle r0..r1 x z0..z1 revolved about the z axis.");
        drag1("r0", app.rv_r0, 0.0625, 0.0, 1e9);
        drag1("r1", app.rv_r1, 0.0625, 0.0, 1e9);
        drag1("z0", app.rv_z0, 0.0625);
        drag1("z1", app.rv_z1, 0.0625);
        ImGui::SliderInt("Quarter turns", &app.rv_quarters, 1, 4);
        static const char* edges[] = {"Outer top", "Outer bottom", "Inner top", "Inner bottom"};
        corner_inputs("rv", app.rv_corner, app.rv_r, app.rv_d2);
        if (app.rv_corner) ImGui::Combo("Edge", &app.rv_edge, edges, 4);
        if (ImGui::Button("Revolve")) {
            const double r0 = app.rv_r0, r1 = app.rv_r1, z0 = app.rv_z0, z1 = app.rv_z1, r = app.rv_r, d2 = app.rv_d2;
            const int q = app.rv_quarters, e = app.rv_edge;
            const auto c = static_cast<wb::Corner>(app.rv_corner);
            app.run("Revolve", [=] { return wb::revolve_rect(r0, r1, z0, z1, q, e, c, r, d2); });
        }
    }

    if (ImGui::CollapsingHeader("Sweep and loft (surfaces)")) {
        ImGui::TextWrapped("Open NURBS surfaces placed at Origin (Primitives): the kernel builds sweeps and lofts as surfaces today, not solids.");
        drag1("Profile radius", app.sw_r, 0.03125, 0.0, 1e9);
        drag3("Path ctrl 1", app.sw_c1);
        drag3("Path ctrl 2", app.sw_c2);
        drag3("Path end", app.sw_end);
        if (ImGui::Button("Sweep circle along path")) {
            const double r = app.sw_r; const Vec3d o = app.origin, a = app.sw_c1, b = app.sw_c2, e = app.sw_end;
            app.run("Sweep", [=] { return wb::sweep_circle(o, r, a, b, e); });
        }
        help("Translational sweep: an exact circle moved along a cubic Bezier path from the origin, without rotating.");
        ImGui::Separator();
        drag1("Radius bottom", app.lf_r0, 0.03125, 0.0, 1e9);
        drag1("Radius middle", app.lf_r1, 0.03125, 0.0, 1e9);
        drag1("Radius top", app.lf_r2, 0.03125, 0.0, 1e9);
        drag1("Height", app.lf_h, 0.0625, 0.0, 1e9);
        if (ImGui::Button("Loft through 3 circles")) {
            const Vec3d o = app.origin; const double a = app.lf_r0, b = app.lf_r1, c = app.lf_r2, h = app.lf_h;
            app.run("Loft", [=] { return wb::loft_circles(o, a, b, c, h); });
        }
    }

    if (ImGui::CollapsingHeader("Freeform fillet")) {
        ImGui::TextWrapped("Block under a freeform top z = h + c x s(y); its edge at x = a rounded by a rolling ball (a procedural blend face).");
        drag1("a", app.ff_a, 0.0625, 0.0, 1e9);
        drag1("L", app.ff_L, 0.0625, 0.0, 1e9);
        drag1("h", app.ff_h, 0.0625, 0.0, 1e9);
        drag1("c (slope)", app.ff_c, 0.03125);
        drag1("Radius", app.ff_r, 0.015625, 0.0, 1e9);
        drag1("End cut (0..1]", app.ff_cut, 0.03125, 0.0, 1.0);
        if (ImGui::Button("Fillet the edge")) {
            const double a = app.ff_a, L = app.ff_L, h = app.ff_h, c = app.ff_c, r = app.ff_r, cut = app.ff_cut;
            app.run("Freeform fillet", [=] { return wb::freeform_fillet(a, L, h, c, r, cut); });
        }
    }

    const Body* A = app.body(app.selA);
    const Body* B = app.body(app.selB);
    if (ImGui::CollapsingHeader("Booleans", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::TextColored(ImVec4(0.45f, 0.7f, 1, 1), "A: %s", A ? A->name.c_str() : "(click a body)");
        ImGui::TextColored(ImVec4(1, 0.65f, 0.3f, 1), "B: %s", B ? B->name.c_str() : "(Shift+click a body)");
        const bool can = A && B && A->solid && B->solid && app.selA != app.selB;
        ImGui::BeginDisabled(!can);
        static const char* labels[3] = {"Union", "Subtract A - B", "Intersect"};
        for (int k = 0; k < 3; ++k) {
            if (k) ImGui::SameLine();
            if (ImGui::Button(labels[k]) && can) {
                const Body a = *A, b = *B;
                app.run(labels[k], [a, b, k] { return wb::boolean(a, b, k); }, {app.selA, app.selB});
            }
        }
        ImGui::EndDisabled();
        ImGui::Checkbox("Hide operands after", &app.hide_operands);
    }

    if (ImGui::CollapsingHeader("Plane cut and move")) {
        ImGui::TextWrapped("On body A. The cut keeps the side n . x <= d.");
        drag3("Normal n", app.cut_n, 0.0625);
        drag1("d", app.cut_d, 0.0625);
        ImGui::BeginDisabled(!(A && A->solid));
        if (ImGui::Button("Cut A by plane") && A) {
            const Body a = *A; const Vec3d n = app.cut_n; const double d = app.cut_d;
            app.run("Plane cut", [=] { return wb::plane_cut(a, n, d); }, {app.selA});
        }
        ImGui::Separator();
        drag3("Translation", app.move_t);
        using SP = ak::op::SignedPermutation;
        auto move = [&](const char* what, SP m) {
            if (!A) return;
            const Body a = *A;
            app.run(what, [=] { return wb::place(a, m, what); }, {app.selA});
        };
        if (ImGui::Button("Translate A")) move("translated", SP::translation(app.move_t));
        if (ImGui::Button("Turn 90 X")) { SP m; m.perm = {0, 2, 1}; m.sign = {1, -1, 1}; move("turned about x", m); }
        ImGui::SameLine();
        if (ImGui::Button("Turn 90 Y")) { SP m; m.perm = {2, 1, 0}; m.sign = {1, 1, -1}; move("turned about y", m); }
        ImGui::SameLine();
        if (ImGui::Button("Turn 90 Z")) move("turned about z", SP::quarter_turn_z());
        if (ImGui::Button("Mirror X")) { SP m; m.sign = {-1, 1, 1}; move("mirrored in x", m); }
        ImGui::SameLine();
        if (ImGui::Button("Mirror Y")) { SP m; m.sign = {1, -1, 1}; move("mirrored in y", m); }
        ImGui::SameLine();
        if (ImGui::Button("Mirror Z")) { SP m; m.sign = {1, 1, -1}; move("mirrored in z", m); }
        help("Exact maps: every control point moves without rounding for dyadic translations (multiples of 1/2^k).");
        ImGui::EndDisabled();
    }

    if (ImGui::CollapsingHeader("View")) {
        ImGui::Checkbox("Edges", &app.show_edges);
        ImGui::SameLine();
        ImGui::Checkbox("Grid", &app.show_grid);
        if (ImGui::Button("Fit view (F)")) app.fit_all();
        ImGui::SameLine();
        if (ImGui::Button("Demo scene")) app.demo();
        ImGui::TextWrapped("Left drag: orbit. Right/middle drag: pan. Wheel: zoom. Click: select A. Shift+click: select B.");
    }
    ImGui::PopItemWidth();
    ImGui::End();
}

void bodies_panel(App& app, float width, float height) {
    const float w = 400;
    ImGui::SetNextWindowPos(ImVec2(width - w, 0), ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(w, height), ImGuiCond_Always);
    ImGui::Begin("Bodies", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);

    if (app.worker.busy()) {
        ImGui::TextColored(ImVec4(1, 0.85f, 0.3f, 1), "Working: %s  %.1f s", app.worker.label().c_str(), app.worker.elapsed());
        ImGui::Separator();
    }

    int remove = -1;
    if (ImGui::BeginChild("list", ImVec2(0, height * 0.3f), ImGuiChildFlags_Borders)) {
        if (app.bodies.empty()) ImGui::TextDisabled("No bodies yet: use the Tools panel.");
        for (std::size_t i = 0; i < app.bodies.size(); ++i) {
            Body& b = app.bodies[i];
            const int idx = static_cast<int>(i);
            ImGui::PushID(idx);
            ImGui::Checkbox("##vis", &b.visible);
            ImGui::SameLine();
            const char* role = idx == app.selA ? "[A] " : idx == app.selB ? "[B] " : "";
            const std::string label = wb::fmt("%s%d  %s%s", role, idx + 1, b.name.c_str(), b.solid ? "" : "  (surface)");
            if (ImGui::Selectable(label.c_str(), idx == app.selA || idx == app.selB, 0, ImVec2(w - 110, 0))) {
                if (ImGui::GetIO().KeyShift) { app.selB = idx; if (app.selA == idx) app.selA = -1; }
                else { app.selA = idx; if (app.selB == idx) app.selB = -1; }
            }
            ImGui::SameLine();
            ImGui::BeginDisabled(app.worker.busy()); // a running job refers to bodies by index
            if (ImGui::SmallButton("Delete")) remove = idx;
            ImGui::EndDisabled();
            ImGui::PopID();
        }
    }
    ImGui::EndChild();
    if (remove >= 0) {
        app.bodies.erase(app.bodies.begin() + remove);
        auto fix = [&](int& s) { s = s == remove ? -1 : s > remove ? s - 1 : s; };
        fix(app.selA);
        fix(app.selB);
    }

    auto inspect = [&](int i, const char* tag, ImVec4 col) {
        const Body* b = app.body(i);
        if (!b) return;
        ImGui::TextColored(col, "%s: %s", tag, b->name.c_str());
        for (const std::string& line : b->info) {
            ImGui::Bullet();
            ImGui::TextWrapped("%s", line.c_str());
        }
        if (b->solid && !b->tess.tris.empty()) {
            ImGui::PushID(tag);
            if (ImGui::SmallButton("Export OBJ")) {
                const std::string path = wb::fmt("akernel_body_%d.obj", i + 1);
                std::ofstream out(path);
                ak::op::write_obj(out, b->tess);
                app.say(static_cast<bool>(out), (out ? "Wrote " : "Could not write ") + path);
            }
            ImGui::PopID();
        }
        ImGui::Spacing();
    };
    if (ImGui::BeginChild("inspector", ImVec2(0, height * 0.38f), ImGuiChildFlags_Borders)) {
        ImGui::PushTextWrapPos(0.0f);
        inspect(app.selA, "A", ImVec4(0.45f, 0.7f, 1, 1));
        inspect(app.selB, "B", ImVec4(1, 0.65f, 0.3f, 1));
        if (app.selA < 0 && app.selB < 0) ImGui::TextDisabled("Select a body to see what the kernel certified.");
        ImGui::PopTextWrapPos();
    }
    ImGui::EndChild();

    ImGui::TextUnformatted("Log");
    if (ImGui::BeginChild("log", ImVec2(0, 0), ImGuiChildFlags_Borders)) {
        ImGui::PushTextWrapPos(0.0f);
        for (const LogLine& l : app.log)
            ImGui::TextColored(l.ok ? ImVec4(0.6f, 0.9f, 0.6f, 1) : ImVec4(1, 0.5f, 0.45f, 1), "%s", l.text.c_str());
        ImGui::PopTextWrapPos();
    }
    ImGui::EndChild();
    ImGui::End();
}

void save_ppm(const char* path, int w, int h) {
    std::vector<unsigned char> px(static_cast<std::size_t>(w) * static_cast<std::size_t>(h) * 3);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, px.data());
    std::ofstream out(path, std::ios::binary);
    out << "P6\n" << w << " " << h << "\n255\n";
    for (int y = h - 1; y >= 0; --y)
        out.write(reinterpret_cast<const char*>(px.data() + static_cast<std::size_t>(y) * static_cast<std::size_t>(w) * 3),
                  static_cast<std::streamsize>(w) * 3);
}

} // namespace

int main(int argc, char** argv) {
    bool demo = false;
    const char* shot = nullptr;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--demo") demo = true;
        else if (a == "--shot" && i + 1 < argc) shot = argv[++i];
    }

    if (!glfwInit()) return 1;
    glfwWindowHint(GLFW_SAMPLES, 4);
    GLFWwindow* window = glfwCreateWindow(1500, 900, "akernel workbench", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        return 1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.IniFilename = nullptr;
    ImGui::StyleColorsDark();
    float xs = 1, ys = 1;
    glfwGetWindowContentScale(window, &xs, &ys);
    if (xs > 1.0f) {
        ImGui::GetStyle().ScaleAllSizes(xs);
        io.FontGlobalScale = xs;
    }
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL2_Init();

    App app;
    if (demo || shot) app.demo();
    app.say(true, "Ready. Every result lists what the kernel certified; refusals show the kernel's reason.");

    double lastx = 0, lasty = 0, pressx = 0, pressy = 0;
    bool dragging = false;
    int frames = 0;
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        app.poll();

        int fw = 0, fh = 0, ww = 0, wh = 0;
        glfwGetFramebufferSize(window, &fw, &fh);
        glfwGetWindowSize(window, &ww, &wh);

        ImGui_ImplOpenGL2_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // Viewport input, unless the UI has the mouse.
        double mx = 0, my = 0;
        glfwGetCursorPos(window, &mx, &my);
        const bool ui = io.WantCaptureMouse;
        const bool l = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        const bool r = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS ||
                       glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS;
        if (!ui && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            dragging = true;
            pressx = mx;
            pressy = my;
        }
        if (dragging && !l) {
            dragging = false;
            if (std::fabs(mx - pressx) + std::fabs(my - pressy) < 4.0) {
                Vec3d o, d;
                app.cam.ray(mx, my, ww, wh, o, d);
                const int hit = app.pick(o, d);
                if (io.KeyShift) { app.selB = hit; if (app.selA == hit) app.selA = -1; }
                else { app.selA = hit; if (app.selB == hit && hit >= 0) app.selB = -1; }
            }
        }
        const double dx = mx - lastx, dy = my - lasty;
        if (!ui || dragging) {
            if (dragging && l) {
                app.cam.yaw -= dx * 0.008;
                app.cam.pitch = std::clamp(app.cam.pitch + dy * 0.008, -1.55, 1.55);
            } else if (r && !ui) {
                Vec3d f, s, u;
                app.cam.frame(f, s, u);
                const double k = 2.0 * app.cam.dist * std::tan(0.5 * app.cam.fovy) / std::max(wh, 1);
                app.cam.target = app.cam.target - s * (dx * k) + u * (dy * k);
            }
            if (!ui && io.MouseWheel != 0.0f) app.cam.dist *= std::pow(0.88, static_cast<double>(io.MouseWheel));
        }
        if (!io.WantCaptureKeyboard && ImGui::IsKeyPressed(ImGuiKey_F)) app.fit_all();
        lastx = mx;
        lasty = my;

        const float W = static_cast<float>(ww), H = static_cast<float>(wh);
        tools_panel(app, H);
        bodies_panel(app, W, H);
        ImGui::Render();

        draw_scene(app, fw, fh);
        ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());

        if (shot && ++frames == 3) {
            save_ppm(shot, fw, fh);
            glfwSetWindowShouldClose(window, 1);
        }
        glfwSwapBuffers(window);
    }

    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
