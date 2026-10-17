// atomsim-viz: the atom, a single fission, and a chain reaction, drawn in real
// time from the same physics library the CLI and tests use.
//
//   atomsim-viz                          interactive
//   atomsim-viz --screenshot 2 out.png [frames] [orbital]
//                                        render a scene, save a PNG, exit
#include "chain.hpp"
#include "constants.hpp"
#include "fission.hpp"
#include "hydrogen.hpp"
#include "nucleus.hpp"

#include "raylib.h"
#include "raymath.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#if defined(__EMSCRIPTEN__)
#include <emscripten/emscripten.h>
#endif

namespace {

// ───────────────────────────── look & feel ──────────────────────────────────

const Color kBackground{9, 11, 20, 255};
const Color kPanel{20, 24, 40, 220};
const Color kText{226, 230, 240, 255};
const Color kDim{140, 148, 170, 255};
const Color kAccent{255, 170, 60, 255};
const Color kProton{235, 80, 80, 255};
const Color kNeutron{150, 170, 200, 255};
const Color kPositive{80, 190, 255, 255};
const Color kNegative{255, 120, 70, 255};

void text(const std::string& s, float x, float y, int size = 20, Color c = kText) {
  DrawText(s.c_str(), int(x), int(y), size, c);
}

std::string fmt(const char* f, double a) {
  char buf[128];
  std::snprintf(buf, sizeof buf, f, a);
  return buf;
}

std::string nuclide(int Z, int A) { return std::to_string(A) + elementSymbol(Z); }

// Perspective camera looking at the origin from `eye`, panned so the origin
// lands at screen x = `screenX` (fraction of the width), clear of the panels.
Camera3D framedCamera(Vector3 eye, float screenX) {
  float dist = Vector3Length(eye);
  float aspect = float(GetScreenWidth()) / float(GetScreenHeight());
  float halfWidth = dist * std::tan(22.5f * DEG2RAD) * aspect;
  Vector3 forward = Vector3Normalize(Vector3Negate(eye));
  Vector3 right = Vector3Normalize(Vector3CrossProduct(forward, {0, 1, 0}));
  Vector3 pan = Vector3Scale(right, -(screenX - 0.5f) * 2 * halfWidth);
  Camera3D cam{};
  cam.position = Vector3Add(eye, pan);
  cam.target = pan;
  cam.up = {0, 1, 0};
  cam.fovy = 45;
  cam.projection = CAMERA_PERSPECTIVE;
  return cam;
}

void panel(float x, float y, float w, float h) {
  DrawRectangleRounded({x, y, w, h}, 0.06f, 8, kPanel);
}

// Points spread through a ball with roughly uniform spacing: random draws
// relaxed by short-range repulsion. Used to draw nucleons.
std::vector<Vector3> packBall(int count, float radius, unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<float> u(-1, 1);
  std::vector<Vector3> p;
  while (int(p.size()) < count) {
    Vector3 v{u(rng), u(rng), u(rng)};
    if (Vector3Length(v) <= 1) p.push_back(Vector3Scale(v, radius));
  }
  const float spacing = 1.8f;
  for (int it = 0; it < 30; ++it)
    for (std::size_t i = 0; i < p.size(); ++i)
      for (std::size_t j = i + 1; j < p.size(); ++j) {
        Vector3 d = Vector3Subtract(p[i], p[j]);
        float len = Vector3Length(d);
        if (len < spacing && len > 1e-4f) {
          Vector3 push = Vector3Scale(d, 0.5f * (spacing - len) / len);
          p[i] = Vector3Add(p[i], push);
          p[j] = Vector3Subtract(p[j], push);
        }
      }
  for (auto& v : p)
    if (Vector3Length(v) > radius) v = Vector3Scale(Vector3Normalize(v), radius);
  return p;
}

// ─────────────────────────────── scene 1: atom ──────────────────────────────

struct Orbital {
  const char* name;
  int n, l;
  RadialState state;
  std::vector<Vector3> points;  // Bohr radii
  std::vector<bool> positive;   // sign of the wavefunction (phase)
};

class AtomScene {
 public:
  AtomScene() {
    for (auto [name, n, l] : {std::tuple{"1s", 1, 0}, {"2s", 2, 0}, {"2p", 2, 1},
                              {"3s", 3, 0}, {"3p", 3, 1}, {"3d", 3, 2}})
      orbitals_.push_back(sample(name, n, l));
  }

  void setOrbital(std::size_t i) { current_ = i % orbitals_.size(); }

  void update(float dt) {
    if (IsKeyPressed(KEY_RIGHT)) current_ = (current_ + 1) % orbitals_.size();
    if (IsKeyPressed(KEY_LEFT)) current_ = (current_ + orbitals_.size() - 1) % orbitals_.size();
    angle_ += dt * 0.25f;
  }

  void draw() {
    const Orbital& o = orbitals_[current_];
    float extent = 3.0f * float(o.n * o.n);
    Camera3D cam = framedCamera(
        {std::cos(angle_) * extent * 2.2f, extent * 0.9f, std::sin(angle_) * extent * 2.2f}, 0.36f);

    BeginBlendMode(BLEND_ADDITIVE);
    for (std::size_t i = 0; i < o.points.size(); ++i) {
      Vector2 s = GetWorldToScreen(o.points[i], cam);
      Color c = o.positive[i] ? kPositive : kNegative;
      DrawRectangleV(s, {2, 2}, Fade(c, 0.35f));
    }
    EndBlendMode();
    Vector2 nucleus = GetWorldToScreen({0, 0, 0}, cam);
    DrawCircleV(nucleus, 3, kProton);

    float w = float(GetScreenWidth());
    panel(w - 380, 80, 360, 300);
    text("Hydrogen atom", w - 360, 95, 30, kAccent);
    text(std::string("orbital ") + o.name + (o.l ? "  (m = 0)" : ""), w - 360, 135);
    text(fmt("E numerical  %.9f eV", o.state.energyEV), w - 360, 170, 20, kText);
    text(fmt("E exact      %.9f eV", schrodingerEnergyEV(o.n, 1.0, phys::kHydrogenReducedMass)),
         w - 360, 195, 20, kDim);
    text(fmt("<r> = %.4f Bohr radii", o.state.meanRadius), w - 360, 225);
    text(fmt("radial nodes: %.0f", o.n - o.l - 1), w - 360, 250);
    text("Each dot: a sampled electron position", w - 360, 285, 10, kDim);
    text("from |psi|^2 (Numerov solution).", w - 360, 300, 10, kDim);
    text("Blue/orange: sign of psi (phase).", w - 360, 315, 10, kDim);
    text("Red dot: the proton, drawn 10,000x", w - 360, 340, 10, kDim);
    text("too big to be visible at all.", w - 360, 355, 10, kDim);
    text("<- ->  change orbital", 20, float(GetScreenHeight()) - 40, 20, kDim);
  }

 private:
  static Orbital sample(const char* name, int n, int l) {
    Orbital o{name, n, l, solveHydrogenic(n, l, 1.0, phys::kHydrogenReducedMass), {}, {}};
    const auto& u = o.state.u;
    std::vector<double> cdf(u.size());
    double sum = 0;
    for (std::size_t i = 0; i < u.size(); ++i) cdf[i] = sum += u[i] * u[i];
    std::mt19937 rng(n * 10 + l);
    std::uniform_real_distribution<double> uni(0, 1);
    // Real m = 0 angular parts: Y00 = 1, Y10 ~ cos, Y20 ~ 3cos^2 - 1.
    auto angular = [l](double c) { return l == 0 ? 1.0 : l == 1 ? c : 1.5 * c * c - 0.5; };
    for (int k = 0; k < 24000; ++k) {
      std::size_t i = std::size_t(std::upper_bound(cdf.begin(), cdf.end(), uni(rng) * sum) - cdf.begin());
      i = std::min(i, u.size() - 1);
      double r = i * o.state.step, cosT;
      do cosT = 2 * uni(rng) - 1;
      while (uni(rng) > angular(cosT) * angular(cosT));
      double phi = 2 * phys::kPi * uni(rng), sinT = std::sqrt(1 - cosT * cosT);
      o.points.push_back({float(r * sinT * std::cos(phi)), float(r * cosT), float(r * sinT * std::sin(phi))});
      o.positive.push_back(u[i] * angular(cosT) >= 0);
    }
    return o;
  }

  std::vector<Orbital> orbitals_;
  std::size_t current_ = 0;
  float angle_ = 0;
};

// ───────────────────────────── scene 2: fission ─────────────────────────────

class FissionScene {
 public:
  FissionScene() : rng_(2026) { restart(); }

  void restart() {
    cn_ = formCompound(92, 235, 0.0253e-6);
    ev_ = simulateFission(cn_, rng_);
    traj_ = coulombTrajectory(ev_, 1e4);
    nucleons_ = packBall(cn_.A, float(nuclearRadius(cn_.A)) * 1.25f, 7);
    // The leftmost A_L nucleons become the light fragment.
    std::vector<int> order(nucleons_.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = int(i);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return nucleons_[a].x < nucleons_[b].x; });
    light_.assign(nucleons_.size(), false);
    isProton_.assign(nucleons_.size(), false);
    for (int k = 0; k < ev_.light.A; ++k) light_[order[k]] = true;
    fragmentIndex_.assign(nucleons_.size(), 0);
    int nl = 0, nh = 0;
    for (std::size_t i = 0; i < nucleons_.size(); ++i) fragmentIndex_[i] = light_[i] ? nl++ : nh++;
    lightHome_ = packBall(ev_.light.A, float(nuclearRadius(ev_.light.A)) * 1.25f, 11);
    heavyHome_ = packBall(ev_.heavy.A, float(nuclearRadius(ev_.heavy.A)) * 1.25f, 13);
    int lp = 0, hp = 0;
    for (std::size_t i = 0; i < nucleons_.size(); ++i) {  // Z protons in each fragment
      if (light_[i]) isProton_[i] = lp++ < ev_.light.Z;
      else isProton_[i] = hp++ < ev_.heavy.Z;
    }
    particles_.clear();
    rings_.clear();
    emitted_ = false;
    clock_ = 0;
  }

  void update(float dt) {
    if (IsKeyPressed(KEY_R) || IsKeyPressed(KEY_SPACE) || clock_ > 14) restart();
    clock_ += dt;
    for (auto& p : particles_) p.pos = Vector3Add(p.pos, Vector3Scale(p.vel, dt));
    for (auto& r : rings_) r.age += dt;
    if (!emitted_ && physicalTime() >= 1e3) emit();
  }

  void draw() {
    const float scission = 2.8f;
    const float D = displaySeparation();
    const double mL = ev_.light.mass, mH = ev_.heavy.mass;
    const float centreL = -D * float(mH / (mL + mH)), centreH = D * float(mL / (mL + mH));
    float zoom = 70 + 0.9f * D;
    Camera3D cam = framedCamera({0, zoom * 0.35f, zoom}, 0.67f);

    BeginMode3D(cam);
    // Incoming neutron.
    if (clock_ < 1.2f) DrawSphere({-60 + 50 * clock_ / 1.2f, 0, 0}, 1.0f, kNeutron);
    for (std::size_t i = 0; i < nucleons_.size(); ++i) {
      Vector3 p = nucleons_[i];
      if (clock_ >= 1.2f) {
        // Before scission: an illustrative stretch from one drop into two.
        // After it, each fragment sits where the computed trajectory puts it.
        float s = std::clamp((clock_ - 1.2f) / (scission - 1.2f), 0.0f, 1.0f);
        float wobble = 0.15f * std::sin(clock_ * 18) * (1 - s);
        Vector3 target = (light_[i] ? lightHome_ : heavyHome_)[fragmentIndex_[i]];
        target.x += light_[i] ? centreL : centreH;
        p = Vector3Lerp(Vector3Scale(p, 1 + wobble), target, s * s * (3 - 2 * s));
      }
      DrawSphereEx(p, 0.85f, 6, 8, isProton_[i] ? kProton : kNeutron);
    }
    for (const auto& p : particles_) {
      DrawSphereEx(p.pos, 1.1f, 6, 8, WHITE);
      DrawLine3D(p.pos, Vector3Subtract(p.pos, Vector3Scale(p.vel, 0.25f)), Fade(WHITE, 0.4f));
    }
    EndMode3D();

    if (clock_ >= scission) {  // label the fragments
      for (auto [x, f] : {std::pair{centreL, &ev_.light}, {centreH, &ev_.heavy}}) {
        float lift = float(nuclearRadius(f->A)) * 1.4f + 3;
        Vector2 s = GetWorldToScreen({x, lift, 0}, cam);
        std::string name = nuclide(f->Z, f->A);
        text(name, s.x - MeasureText(name.c_str(), 20) / 2.0f, s.y - 24, 20, kAccent);
      }
    }
    for (const auto& r : rings_) {  // prompt gamma rays: expanding flashes
      if (r.age < 0 || r.age > 1.2f) continue;
      Vector2 s = GetWorldToScreen(r.at, cam);
      DrawCircleLinesV(s, 20 + 260 * r.age, Fade(YELLOW, 1 - r.age / 1.2f));
    }
    drawPanel();
    drawTimeline();
  }

 private:
  struct Particle { Vector3 pos, vel; };
  struct Ring { Vector3 at; float age; };

  // Display seconds after scission map to physical time on a log scale:
  // 1e-24 s at scission, ten times later every 0.6 s.
  double physicalTime() const {  // zs
    float since = clock_ - 2.8f;
    return since < 0 ? 0 : 1e-3 * std::pow(10.0, since / 0.6);
  }

  TrajectoryPoint stateAt(double t) const {
    if (t <= traj_.front().time) return traj_.front();
    for (std::size_t i = 1; i < traj_.size(); ++i)
      if (t <= traj_[i].time) {
        double a = std::log(std::max(traj_[i - 1].time, 1e-4)), b = std::log(traj_[i].time);
        double f = std::clamp((std::log(t) - a) / (b - a), 0.0, 1.0);
        TrajectoryPoint p = traj_[i - 1];
        auto lerp = [f](double x, double y) { return x + f * (y - x); };
        p.separation = lerp(traj_[i - 1].separation, traj_[i].separation);
        p.kineticLight = lerp(traj_[i - 1].kineticLight, traj_[i].kineticLight);
        p.kineticHeavy = lerp(traj_[i - 1].kineticHeavy, traj_[i].kineticHeavy);
        p.time = t;
        return p;
      }
    return traj_.back();
  }

  // Centre-to-centre distance, drawn on a log scale so 17 fm to 2e5 fm fits.
  float displaySeparation() const {
    double r = clock_ < 2.8f ? ev_.scissionDistance : stateAt(physicalTime()).separation;
    return float(ev_.scissionDistance + 14 * std::log10(r / ev_.scissionDistance));
  }

  void emit() {
    emitted_ = true;
    const float D = displaySeparation();
    const double mL = ev_.light.mass, mH = ev_.heavy.mass;
    std::mt19937 rng(ev_.light.A * 31 + ev_.heavy.Z);
    std::uniform_real_distribution<float> u(-1, 1);
    for (const Fragment* f : {&ev_.light, &ev_.heavy}) {
      float side = f == &ev_.light ? -1.0f : 1.0f;
      Vector3 at{side < 0 ? -D * float(mH / (mL + mH)) : D * float(mL / (mL + mH)), 0, 0};
      for (double e : f->neutronLabEnergies) {
        Vector3 dir{u(rng), u(rng), u(rng)};
        dir = Vector3Normalize(Vector3Add(Vector3Normalize(dir), {side * 0.6f, 0, 0}));
        particles_.push_back({at, Vector3Scale(dir, 22 * float(std::sqrt(e / 2.0)))});
      }
      rings_.push_back({at, -0.3f});
      rings_.push_back({at, -0.6f});
    }
  }

  void drawPanel() {
    float x = 20, y = 80;
    panel(x, y, 430, 470);
    const Fragment &L = ev_.light, &H = ev_.heavy;
    text("n + 235U -> " + nuclide(cn_.Z, cn_.A) + "* -> " + nuclide(L.Z, L.A) + " + " + nuclide(H.Z, H.A),
         x + 15, y + 15, 20, kAccent);
    TrajectoryPoint now = stateAt(physicalTime());
    double t = physicalTime() * 1e-21;
    text(clock_ < 1.2f ? "neutron approaching" : clock_ < 2.8f ? "236U* excited: E* = 6.55 MeV"
                                                               : "t = " + fmt("%.2e s after scission", t),
         x + 15, y + 50);
    text(fmt("separation      %.4g fm", clock_ < 2.8f ? 0.0 : now.separation), x + 15, y + 85, 20, kText);
    text(fmt("KE light        %6.2f MeV", clock_ < 2.8f ? 0.0 : now.kineticLight), x + 15, y + 110);
    text(fmt("KE heavy        %6.2f MeV", clock_ < 2.8f ? 0.0 : now.kineticHeavy), x + 15, y + 135);
    float frac = clock_ < 2.8f ? 0 : float((now.kineticLight + now.kineticHeavy) / ev_.tke);
    DrawRectangle(int(x + 15), int(y + 165), 400, 10, Fade(kDim, 0.3f));
    DrawRectangle(int(x + 15), int(y + 165), int(400 * frac), 10, kAccent);
    text(fmt("Coulomb energy converted to motion: %.1f%%", 100 * frac), x + 15, y + 182, 10, kDim);
    text(fmt("Q-value %.1f MeV", ev_.q) + fmt("   TKE %.1f MeV", ev_.tke), x + 15, y + 205);
    if (emitted_) {
      text(std::to_string(ev_.neutrons()) + " prompt neutrons (lab energies):", x + 15, y + 240);
      std::string list;
      for (const Fragment* f : {&L, &H})
        for (double e : f->neutronLabEnergies) list += fmt("%.2f  ", e);
      text(list + "MeV", x + 15, y + 265, 20, kText);
      text(fmt("prompt gamma rays %.2f MeV", L.gammaEnergy + H.gammaEnergy), x + 15, y + 295);
      text(fmt("fragment spins %.0f", L.spin) + fmt(" and %.0f hbar", H.spin), x + 15, y + 320);
    }
    if (clock_ > 8.5f) {
      text("then beta decay to stability:", x + 15, y + 355);
      float row = y + 380;
      for (const Fragment* f : {&L, &H}) {
        std::string chain = nuclide(f->Z, f->finalA);
        for (const auto& b : f->betaChain) chain += " > " + nuclide(b.Z + 1, f->finalA);
        text(chain, x + 15, row, 20, kAccent);
        row += 25;
      }
    }
    text("Shape change before scission is illustrative; everything after it is computed.",
         x + 15, y + 445, 10, kDim);
  }

  void drawTimeline() {
    float w = float(GetScreenWidth()), h = float(GetScreenHeight());
    text("R  new random fission", 20, h - 40, 20, kDim);
    std::string legend = "red: protons   grey: neutrons   white: emitted neutrons";
    text(legend, w - MeasureText(legend.c_str(), 20) - 20, h - 40, 20, kDim);
  }

  std::mt19937_64 rng_;
  Compound cn_;
  FissionEvent ev_;
  std::vector<TrajectoryPoint> traj_;
  std::vector<Vector3> nucleons_, lightHome_, heavyHome_;
  std::vector<bool> light_, isProton_;
  std::vector<int> fragmentIndex_;
  std::vector<Particle> particles_;
  std::vector<Ring> rings_;
  bool emitted_ = false;
  float clock_ = 0;
};

// ─────────────────────────── scene 3: chain reaction ────────────────────────

class ChainScene {
 public:
  ChainScene() { restart(); }

  void restart() {
    material_ = uranium(kEnrichments[enrichment_]);
    reaction_ = std::make_unique<ChainReaction>(material_, mass_, seed_++, 3000);
    reaction_->inject(40);
    history_.clear();
    flashes_.clear();
  }

  void update(float dt) {
    if (IsKeyPressed(KEY_R)) restart();
    if (IsKeyPressed(KEY_E)) {
      enrichment_ = (enrichment_ + 1) % kEnrichmentCount;
      restart();
    }
    bool changed = false;
    if (IsKeyDown(KEY_UP)) mass_ *= 1 + dt, changed = true;
    if (IsKeyDown(KEY_DOWN)) mass_ /= 1 + dt, changed = true;
    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
      Vector2 m = GetMousePosition();
      Rectangle s = slider();
      if (CheckCollisionPointRec(m, {s.x, s.y - 10, s.width, s.height + 20})) {
        float f = std::clamp((m.x - s.x) / s.width, 0.0f, 1.0f);
        mass_ = kMinMass * std::pow(kMaxMass / kMinMass, f);
        changed = true;
      }
    }
    mass_ = std::clamp(mass_, kMinMass, kMaxMass);
    if (changed) restart();

    // 0.25 ns of physical time per frame, in 5 sub-steps.
    if (reaction_->population() > 0 && reaction_->population() < 1e12)
      for (int k = 0; k < 5; ++k) {
        reaction_->advance(5e-11);
        for (const Flash& f : reaction_->recentFissions())
          if (flashes_.size() < 600) flashes_.push_back({f, 0});
      }
    for (auto& f : flashes_) f.age += dt;
    while (!flashes_.empty() && flashes_.front().age > 0.6f) flashes_.pop_front();
    history_.push_back(std::max(reaction_->population(), 1.0));
    if (history_.size() > 600) history_.erase(history_.begin());
    angle_ += dt * 0.15f;
  }

  void draw() {
    const float scale = 4;  // display units per cm
    float R = float(reaction_->radius()) * scale;
    float dist = std::max(R * 3.4f, 80.0f);
    Camera3D cam = framedCamera({std::cos(angle_) * dist, dist * 0.4f, std::sin(angle_) * dist}, 0.66f);

    BeginMode3D(cam);
    DrawSphereWires({0, 0, 0}, R, 16, 24, Fade(kDim, 0.25f));
    EndMode3D();
    BeginBlendMode(BLEND_ADDITIVE);
    for (const auto& f : flashes_) {
      Vector2 s = GetWorldToScreen({float(f.at.x) * scale, float(f.at.y) * scale, float(f.at.z) * scale}, cam);
      DrawCircleV(s, 4 + 16 * f.age, Fade(kAccent, 0.8f * (1 - f.age / 0.6f)));
    }
    for (const Neutron& n : reaction_->neutrons()) {
      Vector2 s = GetWorldToScreen({float(n.x) * scale, float(n.y) * scale, float(n.z) * scale}, cam);
      // Fast neutrons white-blue, slowed ones violet.
      float hot = std::clamp(float(std::log10(n.energy / 0.05) / std::log10(40.0)), 0.0f, 1.0f);
      Color c{(unsigned char)(150 + 80 * hot), (unsigned char)(110 + 140 * hot), 255, 255};
      DrawRectangleV(s, {3, 3}, Fade(c, 0.85f));
    }
    EndBlendMode();
    drawPanel();
  }

 private:
  static constexpr float kMinMass = 5, kMaxMass = 3000;
  static constexpr double kEnrichments[] = {0.937, 0.5, 0.2, 0.05, 0.0072};
  static constexpr const char* kEnrichmentNames[] = {"93.7% (weapons grade)", "50%", "20% (research reactor)",
                                                     "5% (power reactor fuel)", "0.72% (natural)"};
  static constexpr int kEnrichmentCount = 5;

  struct Glow { Flash at; float age; };

  Rectangle slider() const { return {40, float(GetScreenHeight()) - 70, 420, 8}; }

  void drawPanel() {
    float x = 20, y = 80;
    panel(x, y, 440, 330);
    text("Chain reaction in a bare uranium sphere", x + 15, y + 15, 20, kAccent);
    text(std::string("enrichment ") + kEnrichmentNames[enrichment_], x + 15, y + 50);
    text(fmt("mass %.1f kg", mass_) + fmt("   radius %.2f cm", reaction_->radius()), x + 15, y + 75);
    text(fmt("time %.1f ns", reaction_->time() * 1e9), x + 15, y + 105);
    text(fmt("neutrons %.3g", reaction_->population()), x + 15, y + 130);
    text(fmt("fissions %.3g", reaction_->fissions()), x + 15, y + 155);
    double joules = reaction_->fissions() * 200 * 1.602176634e-13;
    text(fmt("energy released %.3g J", joules), x + 15, y + 180);

    // Growth over the last 20 ns: e-folding time tells sub/supercritical.
    std::string status = "...";
    Color sc = kDim;
    if (history_.size() > 80) {
      double rate = std::log(history_.back() / history_[history_.size() - 80]) / 20e-9;
      if (reaction_->population() <= 1) status = "chain died out: SUBCRITICAL", sc = kPositive;
      else if (rate > 3e6) status = fmt("SUPERCRITICAL: x e every %.0f ns", 1e9 / rate), sc = kNegative;
      else if (rate < -3e6) status = "dying out: SUBCRITICAL", sc = kPositive;
      else status = "near CRITICAL", sc = kAccent;
    }
    text(status, x + 15, y + 215, 20, sc);
    drawHistory(x + 15, y + 245, 410, 70);

    Rectangle s = slider();
    DrawRectangleRec(s, Fade(kDim, 0.4f));
    float f = std::log(mass_ / kMinMass) / std::log(kMaxMass / kMinMass);
    DrawCircleV({s.x + f * s.width, s.y + 4}, 9, kAccent);
    float crit = std::log(52.42f / kMinMass) / std::log(kMaxMass / kMinMass);
    if (enrichment_ == 0) {
      DrawLineEx({s.x + crit * s.width, s.y - 10}, {s.x + crit * s.width, s.y + 18}, 2, kNegative);
      text("Godiva: 52.4 kg", s.x + crit * s.width - 40, s.y - 28, 10, kNegative);
    }
    text("drag / UP DOWN: mass     E: enrichment     R: restart", s.x, s.y + 20, 10, kDim);
    float w = float(GetScreenWidth()), h = float(GetScreenHeight());
    std::string src = "Continuous-energy neutron transport, ENDF/B-VIII.0 cross-sections";
    text(src, w - MeasureText(src.c_str(), 10) - 20, h - 30, 10, kDim);
  }

  void drawHistory(float x, float y, float w, float h) {
    DrawRectangleLines(int(x), int(y), int(w), int(h), Fade(kDim, 0.4f));
    if (history_.size() < 2) return;
    for (std::size_t i = 1; i < history_.size(); ++i) {
      auto py = [&](double v) { return y + h - float(std::log10(v) / 12.0) * h; };
      float x0 = x + w * float(i - 1) / 600, x1 = x + w * float(i) / 600;
      DrawLineV({x0, py(history_[i - 1])}, {x1, py(history_[i])}, kAccent);
    }
    text("neutron population, log scale 1 to 10^12", x + 5, y + 5, 10, kDim);
  }

  Fuel material_ = uranium(0.937);
  std::unique_ptr<ChainReaction> reaction_;
  std::vector<double> history_;
  std::deque<Glow> flashes_;
  float mass_ = 100, angle_ = 0;
  int enrichment_ = 0;
  std::uint64_t seed_ = 1;
};

// ──────────────────────────────── the app ───────────────────────────────────

struct App {
  AtomScene atom;
  FissionScene fission;
  ChainScene chain;
  int scene = 1;
  bool fixedStep = false;  // screenshots: deterministic 1/60 s per frame

  void frame() {
    float dt = fixedStep ? 1 / 60.0f : std::min(GetFrameTime(), 1 / 20.0f);
    if (IsKeyPressed(KEY_ONE)) scene = 0;
    if (IsKeyPressed(KEY_TWO)) scene = 1;
    if (IsKeyPressed(KEY_THREE)) scene = 2;
    for (int i = 0; i < 3; ++i)
      if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && CheckCollisionPointRec(GetMousePosition(), tab(i)))
        scene = i;
    if (scene == 0) atom.update(dt);
    if (scene == 1) fission.update(dt);
    if (scene == 2) chain.update(dt);

    BeginDrawing();
    ClearBackground(kBackground);
    if (scene == 0) atom.draw();
    if (scene == 1) fission.draw();
    if (scene == 2) chain.draw();
    text("atomsim", 20, 20, 30, kText);
    static const char* kNames[] = {"1  Atom", "2  Fission", "3  Chain reaction"};
    for (int i = 0; i < 3; ++i) {
      Rectangle r = tab(i);
      DrawRectangleRounded(r, 0.3f, 6, i == scene ? kAccent : kPanel);
      text(kNames[i], r.x + 14, r.y + 8, 20, i == scene ? kBackground : kText);
    }
    EndDrawing();
  }

  static Rectangle tab(int i) { return {170.0f + i * 225, 18, 210, 36}; }
};

App* gApp = nullptr;
void frameThunk() { gApp->frame(); }

}  // namespace

int main(int argc, char** argv) {
  int shotScene = -1, shotFrames = 0, shotOrbital = 0;
  const char* shotPath = nullptr;
  for (int i = 1; i < argc; ++i)
    if (!std::strcmp(argv[i], "--screenshot") && i + 2 < argc) {
      shotScene = std::atoi(argv[i + 1]) - 1;
      shotPath = argv[i + 2];
      shotFrames = i + 3 < argc ? std::atoi(argv[i + 3]) : 300;
      shotOrbital = i + 4 < argc ? std::atoi(argv[i + 4]) : 0;
    }

  SetConfigFlags(FLAG_MSAA_4X_HINT | FLAG_WINDOW_RESIZABLE);
  InitWindow(1280, 720, "atomsim");
  SetTargetFPS(60);
  App app;
  gApp = &app;

  if (shotScene >= 0) {
    app.scene = shotScene;
    app.fixedStep = true;
    app.atom.setOrbital(std::size_t(shotOrbital));
    for (int f = 0; f < shotFrames; ++f) app.frame();
    TakeScreenshot(shotPath);
    CloseWindow();
    return 0;
  }
#if defined(__EMSCRIPTEN__)
  emscripten_set_main_loop(frameThunk, 0, 1);
#else
  while (!WindowShouldClose()) frameThunk();
#endif
  CloseWindow();
  return 0;
}
