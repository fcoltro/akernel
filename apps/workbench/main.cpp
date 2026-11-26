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
// MCP: the app serves the Model Context Protocol at http://127.0.0.1:7878/mcp
// (loopback only) so an assistant such as Claude Code can work in the same
// scene: create bodies, run the kernel's tools, read results, take
// screenshots (mcp_server.hpp, register_tools below).
//
// Command line: --demo builds a small scene at start; --mcp-port N serves MCP
// on another port; --no-mcp turns it off; --shot file.ppm renders the scene
// once, saves the window as a PPM image and exits (for scripted checks).
#include "kernel_ops.hpp"
#include "mcp_server.hpp"

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl2.h>

#include <GLFW/glfw3.h>

#include <stb_image_write.h>

#include <atomic>
#include <cmath>
#include <cstdlib>
#include <deque>
#include <exception>
#include <fstream>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <numbers>
#include <stdexcept>
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
    int next_id = 1;

    // Work handed to the UI thread by the MCP server thread (see on_ui).
    std::mutex ui_mutex;
    std::vector<std::function<void()>> ui_queue;
    bool closing = false;
    // A screenshot asked for by an MCP client: filled after the next frame is drawn.
    struct Capture {
        int max_width = 900;
        bool full_window = false;
        std::promise<std::vector<unsigned char>> png;
    };
    std::shared_ptr<Capture> capture;
    std::string mcp_status;

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
        if (b.solid) b.color = palette[static_cast<std::size_t>(next_id - 1) % 6];
        b.id = next_id++;
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

    /// Index of the body with this id, or -1.
    int index_of(int id) const {
        for (std::size_t i = 0; i < bodies.size(); ++i)
            if (bodies[i].id == id) return static_cast<int>(i);
        return -1;
    }

    /// Run f on the UI thread (from another thread) and wait for its result.
    template <class F>
    auto on_ui(F f) -> decltype(f()) {
        std::packaged_task<decltype(f())()> task(std::move(f));
        auto fut = task.get_future();
        {
            std::lock_guard<std::mutex> lk(ui_mutex);
            ui_queue.push_back([&task] { task(); });
        }
        return fut.get();
    }
    void drain_ui() {
        std::vector<std::function<void()>> q;
        {
            std::lock_guard<std::mutex> lk(ui_mutex);
            q.swap(ui_queue);
        }
        for (auto& f : q) f();
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

    if (ImGui::CollapsingHeader("View", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Checkbox("Edges", &app.show_edges);
        ImGui::SameLine();
        ImGui::Checkbox("Grid", &app.show_grid);
        if (ImGui::Button("Fit view (F)")) app.fit_all();
        ImGui::SameLine();
        if (ImGui::Button("Demo scene")) app.demo();
        ImGui::TextWrapped("Left drag: orbit. Right/middle drag: pan. Wheel: zoom. Click: select A. Shift+click: select B.");
        if (!app.mcp_status.empty()) {
            ImGui::Separator();
            ImGui::TextWrapped("%s", app.mcp_status.c_str());
            help("Model Context Protocol server for an assistant such as Claude Code, on this computer only. Connect with:\n"
                 "claude mcp add --transport http akernel http://127.0.0.1:7878/mcp");
        }
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
            const std::string label = wb::fmt("%s#%d  %s%s", role, b.id, b.name.c_str(), b.solid ? "" : "  (surface)");
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
                const std::string path = wb::fmt("akernel_body_%d.obj", b->id);
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

// ------------------------------------------------------------ MCP tools
//
// Each tool parses its arguments on the server thread, copies what it needs
// from the app on the UI thread (App::on_ui), runs the kernel on the server
// thread, and hands the new body back to the UI thread, which logs it as
// "Claude: ..." so the person at the window sees what was done.

/// Read pixels of the current frame as a PNG, cropped to the viewport between
/// the panels (or the whole window), downscaled to at most max_width pixels.
std::vector<unsigned char> grab_png(int fw, int fh, int ww, int max_width, bool full_window) {
    const double sc = ww > 0 ? double(fw) / ww : 1.0;
    int x0 = 0, x1 = fw;
    if (!full_window) {
        x0 = std::min(fw - 1, static_cast<int>(330 * sc));
        x1 = std::max(x0 + 1, fw - static_cast<int>(400 * sc));
    }
    const int w = x1 - x0, h = fh;
    std::vector<unsigned char> px(static_cast<std::size_t>(w) * static_cast<std::size_t>(h) * 3);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadBuffer(GL_BACK);
    glReadPixels(x0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, px.data());
    // Box filter down to at most max_width, flipping rows (GL is bottom-up).
    const double f = std::max(1.0, double(w) / std::max(64, max_width));
    const int ow = std::max(1, static_cast<int>(w / f)), oh = std::max(1, static_cast<int>(h / f));
    std::vector<unsigned char> out(static_cast<std::size_t>(ow) * static_cast<std::size_t>(oh) * 3);
    for (int y = 0; y < oh; ++y)
        for (int x = 0; x < ow; ++x) {
            const int sx0 = static_cast<int>(x * f), sx1 = std::max(sx0 + 1, static_cast<int>((x + 1) * f));
            const int sy0 = static_cast<int>(y * f), sy1 = std::max(sy0 + 1, static_cast<int>((y + 1) * f));
            unsigned acc[3] = {0, 0, 0}, n = 0;
            for (int yy = sy0; yy < std::min(sy1, h); ++yy)
                for (int xx = sx0; xx < std::min(sx1, w); ++xx) {
                    const std::size_t i = (static_cast<std::size_t>(h - 1 - yy) * static_cast<std::size_t>(w) + static_cast<std::size_t>(xx)) * 3;
                    for (std::size_t k = 0; k < 3; ++k) acc[k] += px[i + k];
                    ++n;
                }
            const std::size_t o = (static_cast<std::size_t>(y) * static_cast<std::size_t>(ow) + static_cast<std::size_t>(x)) * 3;
            for (std::size_t k = 0; k < 3; ++k) out[o + k] = static_cast<unsigned char>(n ? acc[k] / n : 0);
        }
    std::vector<unsigned char> png;
    stbi_write_png_to_func(
        [](void* ctx, void* data, int size) {
            auto* v = static_cast<std::vector<unsigned char>*>(ctx);
            const auto* p = static_cast<const unsigned char*>(data);
            v->insert(v->end(), p, p + size);
        },
        &png, ow, oh, 3, out.data(), ow * 3);
    return png;
}

using mcp::json;

json schema(json props, std::vector<std::string> required = {}) {
    json s = {{"type", "object"}, {"properties", std::move(props)}};
    if (!required.empty()) s["required"] = required;
    return s;
}
json num(const char* d) { return {{"type", "number"}, {"description", d}}; }
json integer(const char* d) { return {{"type", "integer"}, {"description", d}}; }
json vec3(const char* d) {
    return {{"type", "array"}, {"items", {{"type", "number"}}}, {"minItems", 3}, {"maxItems", 3}, {"description", d}};
}
json choice(std::vector<std::string> v, const char* d) { return {{"type", "string"}, {"enum", v}, {"description", d}}; }

Vec3d v3(const json& a, const char* key, Vec3d def) {
    if (!a.contains(key)) return def;
    const json& j = a[key];
    if (!j.is_array() || j.size() != 3) throw std::runtime_error(std::string(key) + " must be [x, y, z]");
    return {j[0].get<double>(), j[1].get<double>(), j[2].get<double>()};
}
double d1(const json& a, const char* key, double def) { return a.contains(key) ? a[key].get<double>() : def; }
int i1(const json& a, const char* key, int def) { return a.contains(key) ? a[key].get<int>() : def; }
std::string s1(const json& a, const char* key, const char* def) {
    return a.contains(key) ? a[key].get<std::string>() : std::string(def);
}

wb::Corner corner_of(const std::string& c) {
    if (c == "none") return wb::Corner::None;
    if (c == "fillet") return wb::Corner::Fillet;
    if (c == "chamfer") return wb::Corner::Chamfer;
    throw std::runtime_error("corner must be none, fillet or chamfer");
}

std::string describe(const Body& b) {
    std::string s = wb::fmt("#%d %s%s%s\n", b.id, b.name.c_str(), b.solid ? "" : " (open surface)", b.visible ? "" : " [hidden]");
    if (!b.mesh.empty())
        s += wb::fmt("  bounding box [%.6g, %.6g, %.6g] .. [%.6g, %.6g, %.6g]\n", b.mesh.lo[0], b.mesh.lo[1],
                     b.mesh.lo[2], b.mesh.hi[0], b.mesh.hi[1], b.mesh.hi[2]);
    for (const std::string& l : b.info) s += "  " + l + "\n";
    return s;
}

/// Run a kernel call for an MCP client and put its body in the scene.
json kernel_call(App& app, const std::string& label, const std::function<Outcome()>& f, std::vector<int> hide_ids = {}) {
    app.on_ui([&] { app.say(true, "Claude: " + label + " ..."); });
    Outcome o;
    try {
        o = wb::timed(f);
    } catch (const std::exception& e) {
        o = wb::refused(std::string("exception: ") + e.what());
    }
    return app.on_ui([&]() -> json {
        if (!o.ok || !o.body) {
            app.say(false, wb::fmt("Claude: refused: %s  (%.2f s)", o.message.c_str(), o.seconds));
            return mcp::text(wb::fmt("Refused by the kernel (%.2f s): %s", o.seconds, o.message.c_str()), true);
        }
        for (int id : hide_ids) {
            const int i = app.index_of(id);
            if (i >= 0) app.bodies[static_cast<std::size_t>(i)].visible = false;
        }
        const bool first = app.bodies.empty();
        app.add(std::move(*o.body));
        if (first) app.fit_all();
        const Body& b = app.bodies.back();
        app.say(true, wb::fmt("Claude: %s -> #%d  (%.2f s)", o.message.c_str(), b.id, o.seconds));
        return mcp::text(wb::fmt("Created body #%d in %.2f s.\n", b.id, o.seconds) + describe(b));
    });
}

/// A copy of the body with this id (for an operation's operand).
Body body_copy(App& app, int id) {
    return app.on_ui([&] {
        const int i = app.index_of(id);
        if (i < 0) throw std::runtime_error(wb::fmt("no body #%d", id));
        return app.bodies[static_cast<std::size_t>(i)];
    });
}

void register_tools(mcp::Server& srv, App& app) {
    srv.set_instructions(
        "akernel workbench: a live 3D view of the akernel certified B-rep kernel. Bodies are numbered by id "
        "(list_bodies). Units are model units; prefer dyadic values (multiples of 1/2^k such as 0.125) so the "
        "kernel's control points stay exact. Every created body is audited by the kernel (topology, volume, "
        "certified volume enclosure where covered, edge gap bound, watertight mesh); a refused operation returns "
        "the kernel's reason as an error. Results appear live in the person's window; use screenshot to see it.");

    srv.add({"list_bodies", "List the bodies in the scene with their ids, names, visibility and the kernel's audit lines.",
             schema(json::object()), [&app](const json&) {
                 return app.on_ui([&] {
                     std::string s;
                     for (const Body& b : app.bodies) s += describe(b);
                     if (s.empty()) s = "The scene is empty.";
                     const Body* a = app.body(app.selA);
                     const Body* b = app.body(app.selB);
                     s += wb::fmt("Selection in the window: A = %s, B = %s\n", a ? ("#" + std::to_string(a->id)).c_str() : "none",
                                  b ? ("#" + std::to_string(b->id)).c_str() : "none");
                     return mcp::text(s);
                 });
             }});

    srv.add({"create_primitive",
             "Create a solid primitive. box: origin + size. bump_box: box whose top has a quadratic hill of height "
             "'bump'. spline_box: box whose top is a quadratic B-spline wave with spans_u x spans_v spans and "
             "amplitude 'wave'. rounded_box: box with all 12 edges rounded by 'radius' (exact cylinders and sphere "
             "octants).",
             schema({{"kind", choice({"box", "bump_box", "spline_box", "rounded_box"}, "primitive")},
                     {"origin", vec3("lowest corner, default [0,0,0]")},
                     {"size", vec3("extent in x, y, z, default [1,1,0.5]")},
                     {"bump", num("bump_box hill height, default 0.25")},
                     {"radius", num("rounded_box radius, default 0.125")},
                     {"spans_u", integer("spline_box spans in x, default 2")},
                     {"spans_v", integer("spline_box spans in y, default 2")},
                     {"wave", num("spline_box wave amplitude, default 0.0625")}},
                    {"kind"}),
             [&app](const json& a) {
                 const std::string k = a.at("kind").get<std::string>();
                 const Vec3d o = v3(a, "origin", {0, 0, 0}), s = v3(a, "size", {1, 1, 0.5});
                 std::function<Outcome()> f;
                 if (k == "box") f = [=] { return wb::box(o, s); };
                 else if (k == "bump_box") { const double b = d1(a, "bump", 0.25); f = [=] { return wb::bump_box(o, s, b); }; }
                 else if (k == "spline_box") {
                     const int m = i1(a, "spans_u", 2), n = i1(a, "spans_v", 2);
                     const double w = d1(a, "wave", 0.0625);
                     f = [=] { return wb::spline_box(o, s, m, n, w); };
                 } else if (k == "rounded_box") { const double r = d1(a, "radius", 0.125); f = [=] { return wb::rounded_box(o, s, r); }; }
                 else throw std::runtime_error("unknown kind " + k);
                 return kernel_call(app, k, f);
             }});

    srv.add({"extrude",
             "Extrude the rectangle width x height (in x, y at origin) by length along z. corner: none, fillet "
             "(radius) or chamfer (setbacks radius on the x = width face and setback_b on the y = height face) on "
             "the edge at (width, height). With step_a1 and step_b1 the profile is an L-step instead (without the "
             "notch x > step_a1, y > step_b1) and its concave inside edge gets the fillet or chamfer.",
             schema({{"origin", vec3("default [0,0,0]")}, {"width", num("default 1")}, {"height", num("default 0.5")},
                     {"length", num("default 1")}, {"corner", choice({"none", "fillet", "chamfer"}, "default none")},
                     {"radius", num("fillet radius or first chamfer setback, default 0.125")},
                     {"setback_b", num("second chamfer setback, default = radius")},
                     {"step_a1", num("L-step notch x")}, {"step_b1", num("L-step notch y")}}),
             [&app](const json& a) {
                 const Vec3d o = v3(a, "origin", {0, 0, 0});
                 const double W = d1(a, "width", 1), H = d1(a, "height", 0.5), L = d1(a, "length", 1);
                 const wb::Corner c = corner_of(s1(a, "corner", "none"));
                 const double r = d1(a, "radius", 0.125), d2 = d1(a, "setback_b", r);
                 if (a.contains("step_a1") || a.contains("step_b1")) {
                     const double a1 = d1(a, "step_a1", W / 2), b1 = d1(a, "step_b1", H / 2);
                     return kernel_call(app, "extrude step", [=] { return wb::extrude_step(o, W, H, L, a1, b1, c, r, d2); });
                 }
                 return kernel_call(app, "extrude", [=] { return wb::extrude_rect(o, W, H, L, c, r, d2); });
             }});

    srv.add({"revolve",
             "Revolve the rectangle r0..r1 x z0..z1 about the z axis by quarter_turns (1-4; 4 is a closed ring). "
             "Optionally fillet or chamfer one of its circular edges.",
             schema({{"r0", num("inner radius > 0, default 0.25")}, {"r1", num("outer radius, default 0.5")},
                     {"z0", num("default 0")}, {"z1", num("default 0.5")}, {"quarter_turns", integer("1..4, default 4")},
                     {"corner", choice({"none", "fillet", "chamfer"}, "default none")},
                     {"edge", choice({"outer_top", "outer_bottom", "inner_top", "inner_bottom"}, "default outer_top")},
                     {"radius", num("fillet radius or chamfer setback on the cylinder, default 0.125")},
                     {"setback_b", num("chamfer setback on the flat face, default = radius")}}),
             [&app](const json& a) {
                 const double r0 = d1(a, "r0", 0.25), r1 = d1(a, "r1", 0.5), z0 = d1(a, "z0", 0), z1 = d1(a, "z1", 0.5);
                 const int q = i1(a, "quarter_turns", 4);
                 if (q < 1 || q > 4) throw std::runtime_error("quarter_turns must be 1..4");
                 const wb::Corner c = corner_of(s1(a, "corner", "none"));
                 const std::string e = s1(a, "edge", "outer_top");
                 const int ei = e == "outer_top" ? 0 : e == "outer_bottom" ? 1 : e == "inner_top" ? 2 : e == "inner_bottom" ? 3 : -1;
                 if (ei < 0) throw std::runtime_error("bad edge");
                 const double r = d1(a, "radius", 0.125), d2 = d1(a, "setback_b", r);
                 return kernel_call(app, "revolve", [=] { return wb::revolve_rect(r0, r1, z0, z1, q, ei, c, r, d2); });
             }});

    srv.add({"sweep_surface",
             "Open NURBS surface: an exact circle of 'radius' at origin moved without rotation along a cubic Bezier "
             "path from origin through origin+ctrl1, origin+ctrl2 to origin+end (translational sweep).",
             schema({{"origin", vec3("default [0,0,0]")}, {"radius", num("default 0.25")}, {"ctrl1", vec3("default [0,0.5,0.5]")},
                     {"ctrl2", vec3("default [0.5,0.5,1]")}, {"end", vec3("default [1,0,1.5]")}}),
             [&app](const json& a) {
                 const Vec3d o = v3(a, "origin", {0, 0, 0}), c1 = v3(a, "ctrl1", {0, 0.5, 0.5}), c2 = v3(a, "ctrl2", {0.5, 0.5, 1}),
                             e = v3(a, "end", {1, 0, 1.5});
                 const double r = d1(a, "radius", 0.25);
                 return kernel_call(app, "sweep", [=] { return wb::sweep_circle(o, r, c1, c2, e); });
             }});

    srv.add({"loft_surface", "Open NURBS surface skinned through three horizontal circles about the vertical axis through origin.",
             schema({{"origin", vec3("default [0,0,0]")}, {"r_bottom", num("default 0.5")}, {"r_middle", num("default 0.25")},
                     {"r_top", num("default 0.375")}, {"height", num("default 1")}}),
             [&app](const json& a) {
                 const Vec3d o = v3(a, "origin", {0, 0, 0});
                 const double r0 = d1(a, "r_bottom", 0.5), r1 = d1(a, "r_middle", 0.25), r2 = d1(a, "r_top", 0.375), h = d1(a, "height", 1);
                 return kernel_call(app, "loft", [=] { return wb::loft_circles(o, r0, r1, r2, h); });
             }});

    srv.add({"freeform_fillet",
             "A block under the freeform top z = h + slope*x*s(y) over [0,a] x [0,L] (s the smoothstep), with its "
             "edge at x = a rounded by a rolling ball of 'radius' (a procedural blend face with a certified spine). "
             "end_cut in (0,1] ends the block at y = end_cut*L.",
             schema({{"a", num("default 1")}, {"L", num("default 2")}, {"h", num("default 0.5")}, {"slope", num("default 0.25")},
                     {"radius", num("default 0.125")}, {"end_cut", num("default 1")}}),
             [&app](const json& a) {
                 const double A = d1(a, "a", 1), L = d1(a, "L", 2), h = d1(a, "h", 0.5), c = d1(a, "slope", 0.25),
                              r = d1(a, "radius", 0.125), cut = d1(a, "end_cut", 1);
                 return kernel_call(app, "freeform fillet", [=] { return wb::freeform_fillet(A, L, h, c, r, cut); });
             }});

    srv.add({"boolean",
             "Boolean of two solids by id: union, difference (a minus b) or intersection. Uses the procedural-face "
             "Booleans when an operand has blend faces. Operands are hidden afterwards unless hide_operands is false.",
             schema({{"a", integer("id of the first solid")}, {"b", integer("id of the second solid")},
                     {"op", choice({"union", "difference", "intersection"}, "operation")},
                     {"hide_operands", {{"type", "boolean"}, {"description", "default true"}}}},
                    {"a", "b", "op"}),
             [&app](const json& a) {
                 const int ia = a.at("a").get<int>(), ib = a.at("b").get<int>();
                 const std::string op = a.at("op").get<std::string>();
                 const int k = op == "union" ? 0 : op == "difference" ? 1 : op == "intersection" ? 2 : -1;
                 if (k < 0) throw std::runtime_error("op must be union, difference or intersection");
                 Body A = body_copy(app, ia), B = body_copy(app, ib);
                 std::vector<int> hide;
                 if (a.value("hide_operands", true)) hide = {ia, ib};
                 return kernel_call(app, wb::fmt("%s of #%d and #%d", op.c_str(), ia, ib), [&] { return wb::boolean(A, B, k); }, hide);
             }});

    srv.add({"plane_cut", "Cut a solid by the plane normal . x = d, keeping the side normal . x <= d (certified section curves).",
             schema({{"id", integer("solid id")}, {"normal", vec3("plane normal")}, {"d", num("plane offset")},
                     {"hide_operand", {{"type", "boolean"}, {"description", "default true"}}}},
                    {"id", "normal", "d"}),
             [&app](const json& a) {
                 const int id = a.at("id").get<int>();
                 const Vec3d n = v3(a, "normal", {0, 0, 1});
                 const double d = a.at("d").get<double>();
                 Body A = body_copy(app, id);
                 std::vector<int> hide;
                 if (a.value("hide_operand", true)) hide = {id};
                 return kernel_call(app, wb::fmt("plane cut of #%d", id), [&] { return wb::plane_cut(A, n, d); }, hide);
             }});

    srv.add({"transform",
             "Move a solid exactly (a new body): translate by 'vector', quarter turns about an axis through the "
             "origin, or mirror in a coordinate plane through the origin.",
             schema({{"id", integer("solid id")},
                     {"op", choice({"translate", "turn_x", "turn_y", "turn_z", "mirror_x", "mirror_y", "mirror_z"}, "map")},
                     {"vector", vec3("translation, for translate")},
                     {"hide_operand", {{"type", "boolean"}, {"description", "default true"}}}},
                    {"id", "op"}),
             [&app](const json& a) {
                 const int id = a.at("id").get<int>();
                 const std::string op = a.at("op").get<std::string>();
                 using SP = ak::op::SignedPermutation;
                 SP m;
                 if (op == "translate") m = SP::translation(v3(a, "vector", {0, 0, 0}));
                 else if (op == "turn_x") { m.perm = {0, 2, 1}; m.sign = {1, -1, 1}; }
                 else if (op == "turn_y") { m.perm = {2, 1, 0}; m.sign = {1, 1, -1}; }
                 else if (op == "turn_z") m = SP::quarter_turn_z();
                 else if (op == "mirror_x") m.sign = {-1, 1, 1};
                 else if (op == "mirror_y") m.sign = {1, -1, 1};
                 else if (op == "mirror_z") m.sign = {1, 1, -1};
                 else throw std::runtime_error("unknown op " + op);
                 Body A = body_copy(app, id);
                 std::vector<int> hide;
                 if (a.value("hide_operand", true)) hide = {id};
                 return kernel_call(app, op + wb::fmt(" of #%d", id), [&] { return wb::place(A, m, ""); }, hide);
             }});

    srv.add({"edit_scene",
             "Change the scene: show, hide or delete bodies by id, select bodies as A and B in the window (as the "
             "person would by clicking), clear everything, or load the demo scene.",
             schema({{"action", choice({"show", "hide", "delete", "select", "clear", "demo"}, "what to do")},
                     {"ids", {{"type", "array"}, {"items", {{"type", "integer"}}}, {"description", "body ids (select: [A] or [A, B])"}}}},
                    {"action"}),
             [&app](const json& a) {
                 const std::string act = a.at("action").get<std::string>();
                 std::vector<int> ids;
                 if (a.contains("ids")) ids = a["ids"].get<std::vector<int>>();
                 return app.on_ui([&]() -> json {
                     if (act == "demo") {
                         app.demo();
                         app.say(true, "Claude: loaded the demo scene");
                         return mcp::text("Demo scene loaded.");
                     }
                     if ((act == "delete" || act == "clear") && app.worker.busy())
                         return mcp::text("The window is running an operation; try again when it finishes.", true);
                     if (act == "clear") {
                         app.bodies.clear();
                         app.selA = app.selB = -1;
                         app.say(true, "Claude: cleared the scene");
                         return mcp::text("Scene cleared.");
                     }
                     std::string missing;
                     for (int id : ids)
                         if (app.index_of(id) < 0) missing += " #" + std::to_string(id);
                     if (!missing.empty()) return mcp::text("No such bodies:" + missing, true);
                     if (act == "select") {
                         app.selA = ids.size() > 0 ? app.index_of(ids[0]) : -1;
                         app.selB = ids.size() > 1 ? app.index_of(ids[1]) : -1;
                     } else if (act == "show" || act == "hide") {
                         for (int id : ids) app.bodies[static_cast<std::size_t>(app.index_of(id))].visible = act == "show";
                     } else if (act == "delete") {
                         for (int id : ids) {
                             const int i = app.index_of(id);
                             app.bodies.erase(app.bodies.begin() + i);
                             auto fix = [&](int& s) { s = s == i ? -1 : s > i ? s - 1 : s; };
                             fix(app.selA);
                             fix(app.selB);
                         }
                     } else {
                         return mcp::text("unknown action " + act, true);
                     }
                     app.say(true, "Claude: " + act + " " + a.value("ids", json::array()).dump());
                     return mcp::text("Done.");
                 });
             }});

    srv.add({"view",
             "Set the camera: a preset direction (iso, top, bottom, front, back, left, right), fit to the visible "
             "bodies, or explicit yaw/pitch (radians, z up), distance and target. Also toggles edges and grid.",
             schema({{"preset", choice({"iso", "top", "bottom", "front", "back", "left", "right"}, "view direction")},
                     {"fit", {{"type", "boolean"}, {"description", "fit the visible bodies (default true with a preset)"}}},
                     {"yaw", num("radians")}, {"pitch", num("radians, -1.55..1.55")}, {"distance", num("from target")},
                     {"target", vec3("look-at point")}, {"edges", {{"type", "boolean"}}}, {"grid", {{"type", "boolean"}}}}),
             [&app](const json& a) {
                 return app.on_ui([&] {
                     Camera& c = app.cam;
                     if (a.contains("preset")) {
                         const std::string p = a["preset"].get<std::string>();
                         const double pi = std::numbers::pi;
                         if (p == "iso") { c.yaw = -pi / 4; c.pitch = 0.6155; }
                         else if (p == "top") { c.yaw = -pi / 2; c.pitch = 1.55; }
                         else if (p == "bottom") { c.yaw = -pi / 2; c.pitch = -1.55; }
                         else if (p == "front") { c.yaw = -pi / 2; c.pitch = 0; }
                         else if (p == "back") { c.yaw = pi / 2; c.pitch = 0; }
                         else if (p == "right") { c.yaw = 0; c.pitch = 0; }
                         else if (p == "left") { c.yaw = pi; c.pitch = 0; }
                     }
                     if (a.contains("yaw")) c.yaw = a["yaw"].get<double>();
                     if (a.contains("pitch")) c.pitch = std::clamp(a["pitch"].get<double>(), -1.55, 1.55);
                     if (a.value("fit", a.contains("preset"))) app.fit_all();
                     if (a.contains("distance")) c.dist = std::max(1e-6, a["distance"].get<double>());
                     if (a.contains("target")) c.target = v3(a, "target", c.target);
                     if (a.contains("edges")) app.show_edges = a["edges"].get<bool>();
                     if (a.contains("grid")) app.show_grid = a["grid"].get<bool>();
                     return mcp::text(wb::fmt("Camera: yaw %.4g, pitch %.4g, distance %.4g, target [%.4g, %.4g, %.4g]", c.yaw,
                                              c.pitch, c.dist, c.target.x, c.target.y, c.target.z));
                 });
             }});

    srv.add({"screenshot", "A PNG of the 3D viewport as the person sees it now (or the whole window with the panels).",
             schema({{"max_width", integer("downscale to at most this many pixels wide, default 900")},
                     {"full_window", {{"type", "boolean"}, {"description", "include the tool panels, default false"}}}}),
             [&app](const json& a) {
                 auto cap = std::make_shared<App::Capture>();
                 cap->max_width = std::clamp(i1(a, "max_width", 900), 64, 4096);
                 cap->full_window = a.value("full_window", false);
                 auto fut = cap->png.get_future();
                 app.on_ui([&] {
                     if (app.closing) cap->png.set_value({});
                     else app.capture = cap;
                 });
                 if (fut.wait_for(std::chrono::seconds(10)) != std::future_status::ready)
                     return mcp::text("The window did not draw a frame in 10 s (is it minimised?)", true);
                 const std::vector<unsigned char> png = fut.get();
                 if (png.empty()) return mcp::text("No image (the window is closing).", true);
                 return mcp::image_png(png, "");
             }});

    srv.add({"export_obj", "Write the kernel's watertight mesh of a solid as a Wavefront OBJ file (one group per face, edges as polylines).",
             schema({{"id", integer("solid id")}, {"path", {{"type", "string"}, {"description", "file path on the person's computer"}}}},
                    {"id", "path"}),
             [&app](const json& a) {
                 const int id = a.at("id").get<int>();
                 const std::string path = a.at("path").get<std::string>();
                 const Body b = body_copy(app, id);
                 if (b.tess.tris.empty()) return mcp::text("This body has no kernel mesh (an open surface, or tessellation was refused).", true);
                 std::ofstream out(path);
                 ak::op::write_obj(out, b.tess);
                 if (!out) return mcp::text("Could not write " + path, true);
                 app.on_ui([&] { app.say(true, "Claude: wrote " + path); });
                 return mcp::text(wb::fmt("Wrote #%d to %s (%d triangles).", id, path.c_str(), static_cast<int>(b.tess.tris.size())));
             }});

    srv.add({"get_log", "The workbench log, newest first: results and refusals of everything done in the window or by MCP.",
             schema({{"limit", integer("lines, default 30")}}), [&app](const json& a) {
                 const int n = std::max(1, i1(a, "limit", 30));
                 return app.on_ui([&] {
                     std::string s;
                     int k = 0;
                     for (const LogLine& l : app.log) {
                         if (k++ >= n) break;
                         s += (l.ok ? "ok   " : "FAIL ") + l.text + "\n";
                     }
                     return mcp::text(s.empty() ? "The log is empty." : s);
                 });
             }});
}

} // namespace

int main(int argc, char** argv) {
    bool demo = false, serve_mcp = true;
    int mcp_port = 7878;
    const char* shot = nullptr;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--demo") demo = true;
        else if (a == "--shot" && i + 1 < argc) shot = argv[++i];
        else if (a == "--no-mcp") serve_mcp = false;
        else if (a == "--mcp-port" && i + 1 < argc) mcp_port = std::atoi(argv[++i]);
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

    mcp::Server server;
    register_tools(server, app);
    if (serve_mcp && !shot) {
        std::string err;
        if (server.start(mcp_port, err)) {
            app.mcp_status = wb::fmt("MCP: http://127.0.0.1:%d/mcp", mcp_port);
            app.say(true, app.mcp_status + " (Claude Code: claude mcp add --transport http akernel " +
                              wb::fmt("http://127.0.0.1:%d/mcp)", mcp_port));
        } else {
            app.mcp_status = "MCP off: " + err;
            app.say(false, app.mcp_status);
        }
    }

    double lastx = 0, lasty = 0, pressx = 0, pressy = 0;
    bool dragging = false;
    int frames = 0;
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        app.poll();
        app.drain_ui();

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
        if (app.capture && !app.capture->full_window) {
            app.capture->png.set_value(grab_png(fw, fh, ww, app.capture->max_width, false));
            app.capture.reset();
        }
        ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
        if (app.capture) {
            app.capture->png.set_value(grab_png(fw, fh, ww, app.capture->max_width, true));
            app.capture.reset();
        }

        if (shot && ++frames == 3) {
            save_ppm(shot, fw, fh);
            glfwSetWindowShouldClose(window, 1);
        }
        glfwSwapBuffers(window);
    }

    // Stop the MCP server; a request in flight may still be waiting for the UI thread.
    app.closing = true;
    if (app.capture) {
        app.capture->png.set_value({});
        app.capture.reset();
    }
    std::atomic<bool> stopped{false};
    std::thread stopper([&] {
        server.stop();
        stopped = true;
    });
    while (!stopped) {
        app.drain_ui();
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    stopper.join();

    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
