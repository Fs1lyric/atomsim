// atomsim: an atom, and the splitting of its nucleus, from measured data and
// established theory. See README.md for the models and their limits.
#include "constants.hpp"
#include "electrons.hpp"
#include "fission.hpp"
#include "hydrogen.hpp"
#include "nucleus.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <random>
#include <string>
#include <thread>
#include <vector>

namespace {

struct Options {
  int targetZ = 92, targetA = 235;
  double neutronEnergy = 0.0253e-6;  // thermal: 0.0253 eV
  int events = 100000;
  unsigned long long seed = 2026;
  int threads = 0;      // 0 = all cores
  std::string csvDir;   // write plot data here when set
};

void usage() {
  std::puts(
      "usage: atomsim [--target U235|U238|Pu239] [--energy MeV] [--events N] [--seed N]\n"
      "               [--threads N] [--csv DIR]\n"
      "  --energy   incoming neutron kinetic energy in MeV (default 2.53e-8, thermal)\n"
      "  --events   Monte Carlo fission events (default 100000)\n"
      "  --threads  worker threads for the Monte Carlo (default: all cores)\n"
      "  --csv      also write the plotted data as CSV files into DIR");
}

Options parse(int argc, char** argv) {
  Options o;
  for (int i = 1; i < argc; ++i) {
    auto next = [&]() -> const char* {
      if (i + 1 >= argc) {
        usage();
        std::exit(2);
      }
      return argv[++i];
    };
    if (!std::strcmp(argv[i], "--target")) {
      std::string t = next();
      if (t == "U235") o.targetZ = 92, o.targetA = 235;
      else if (t == "U238") o.targetZ = 92, o.targetA = 238;
      else if (t == "Pu239") o.targetZ = 94, o.targetA = 239;
      else { usage(); std::exit(2); }
    } else if (!std::strcmp(argv[i], "--energy")) {
      o.neutronEnergy = std::atof(next());
    } else if (!std::strcmp(argv[i], "--events")) {
      o.events = std::max(1, std::atoi(next()));
    } else if (!std::strcmp(argv[i], "--threads")) {
      o.threads = std::max(0, std::atoi(next()));
    } else if (!std::strcmp(argv[i], "--csv")) {
      o.csvDir = next();
    } else if (!std::strcmp(argv[i], "--seed")) {
      o.seed = std::strtoull(next(), nullptr, 10);
    } else {
      usage();
      std::exit(argv[i] == std::string("--help") ? 0 : 2);
    }
  }
  return o;
}

std::string nuclide(int Z, int A) { return std::to_string(A) + elementSymbol(Z); }

// A horizontal bar exactly `width` terminal columns wide.
std::string bar(double fraction, int width) {
  static const char* kEighths[] = {"", "▏", "▎", "▍", "▌", "▋", "▊", "▉"};
  fraction = std::fmax(0.0, std::fmin(1.0, fraction));
  int eighths = int(std::lround(fraction * width * 8));
  int full = eighths / 8, part = eighths % 8;
  std::string s;
  for (int i = 0; i < full; ++i) s += "█";
  s += kEighths[part];
  return s + std::string(width - full - (part ? 1 : 0), ' ');
}

void heading(const char* title) {
  std::printf("\n\033[1m══ %s ", title);
  for (int i = int(std::strlen(title)); i < 72; ++i) std::printf("═");
  std::printf("\033[0m\n");
}

void section(const char* title) { std::printf("\n\033[1m── %s\033[0m\n", title); }

// ───────────────────────────── Part 1: the atom ─────────────────────────────

void hydrogenAtom() {
  section("1.1  Hydrogen: radial Schrödinger equation solved numerically (Numerov shooting)");
  const double mu = phys::kHydrogenReducedMass;
  std::printf("  state   E numerical (eV)   E exact (eV)      error (eV)    <r> (a0)\n");
  const char letters[] = "spd";
  for (int n = 1; n <= 4; ++n)
    for (int l = 0; l < n && l <= 2; ++l) {
      RadialState s = solveHydrogenic(n, l, 1.0, mu);
      double exact = schrodingerEnergyEV(n, 1.0, mu);
      std::printf("  %d%c   %16.9f  %14.9f  %12.2e  %9.4f\n", n, letters[l], s.energyEV, exact,
                  s.energyEV - exact, s.meanRadius);
    }

  std::printf("\n  Radial probability r²|R(r)|² for r = 0..25 Bohr radii:\n");
  static const char* kLevels[] = {" ", "▁", "▂", "▃", "▄", "▅", "▆", "▇", "█"};
  for (auto [n, l] : {std::pair{1, 0}, {2, 0}, {2, 1}, {3, 0}, {3, 2}}) {
    RadialState s = solveHydrogenic(n, l, 1.0, mu);
    const int columns = 64;
    double peak = 0, density[columns];
    for (int c = 0; c < columns; ++c) {
      size_t i = size_t((c + 0.5) * 25.0 / columns / s.step);
      density[c] = i < s.u.size() ? s.u[i] * s.u[i] : 0.0;
      peak = std::fmax(peak, density[c]);
    }
    std::printf("  %d%c │", n, letters[l]);
    for (double d : density) std::printf("%s", kLevels[int(std::lround(8 * d / peak))]);
    std::printf("│\n");
  }

  const double muMeV = phys::kHydrogenReducedMassMeV;
  double fs = diracEnergyEV(2, 1.5, 1.0, muMeV) - diracEnergyEV(2, 0.5, 1.0, muMeV);
  std::printf(
      "\n  Relativistic (Dirac) correction: 1s½ %.6f eV, 2p³⁄₂ - 2p½ splitting %.4f GHz\n"
      "  (measured 10.969 GHz; the rest, and the 1.058 GHz Lamb shift, is QED — not modelled)\n",
      diracEnergyEV(1, 0.5, 1.0, muMeV), fs / phys::kPlanckEVs / 1e9);
}

void heavyAtom(int Z, int A) {
  char title[96];
  std::snprintf(title, sizeof title, "1.2  The %s atom", nuclide(Z, A).c_str());
  section(title);
  double be = bindingEnergy(Z, A), radius = nuclearRadius(A);
  std::printf("  Nucleus     %d protons + %d neutrons, radius R = 1.2·A^⅓ = %.2f fm\n", Z, A - Z,
              radius);
  std::printf("  Binding     %.2f MeV (%s), %.3f MeV per nucleon\n", be,
              isMeasured(Z, A) ? "AME2020" : "model", be / A);
  std::printf("              textbook SEMF %.2f MeV, liquid drop + shells %.2f MeV\n",
              semfBindingEnergy(Z, A), shellBindingEnergy(Z, A));
  std::printf("  Electrons   %d in %s\n", Z, electronConfiguration(Z).c_str());
  std::printf("  Scale       atomic radius ≈ 156 pm = 1.56×10⁵ fm, ~%.0f× the nuclear radius\n",
              1.56e5 / radius);
  std::printf("  Coulomb     repulsion between all %d protons: %.0f MeV, held by the strong force\n",
              Z, 0.6 * phys::kCoulomb * Z * (Z - 1) / radius);

  std::printf("\n  Why heavy nuclei can release energy by splitting — binding per nucleon:\n");
  for (auto [z, a] : {std::pair{2, 4}, {6, 12}, {26, 56}, {82, 208}, {92, 235}}) {
    double perA = bindingEnergy(z, a) / a;
    std::printf("  %-6s %6.3f MeV %s\n", nuclide(z, a).c_str(), perA,
                bar((perA - 6.5) / 2.4, 40).c_str());
  }
  std::printf("  Peak near iron. Two mid-mass fragments bind ~0.9 MeV/nucleon more than U.\n");
}

// ────────────────────────── Part 2: splitting it ────────────────────────────

void capture(const Compound& cn) {
  section("2.1  Neutron capture and the fission barrier");
  std::printf("  n + %s → %s*   (E_n = %.4g MeV)\n", nuclide(cn.targetZ, cn.targetA).c_str(),
              nuclide(cn.Z, cn.A).c_str(), cn.neutronEnergy);
  std::printf("  Neutron separation energy S_n  %7.3f MeV  (AME2020 masses)\n",
              cn.separationEnergy);
  std::printf("  Excitation energy E*           %7.3f MeV  = S_n + E_n·A/(A+1)\n", cn.excitation);
  std::printf("  Fissility x = E_C / 2E_S       %7.3f\n", cn.fissility);
  std::printf("  Liquid-drop fission barrier    %7.3f MeV  (Cohen-Swiatecki / Myers-Swiatecki)\n",
              cn.barrier);
  std::printf("  Hill-Wheeler transmission      %7.3g\n", cn.transmission);

  std::printf("\n  Transmission vs neutron energy:\n");
  for (double en : {0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0}) {
    Compound c = formCompound(cn.targetZ, cn.targetA, en + 0.0253e-6);
    std::printf("  %4.2f MeV │%s│ %.3g\n", en, bar(c.transmission, 40).c_str(),
                c.transmission);
  }
}

void filmStrip(const FissionEvent& ev, const std::vector<TrajectoryPoint>& traj) {
  std::printf("\n  The fragments flying apart (log scale either side of the centre of mass):\n");
  const int half = 34;
  double mL = ev.light.mass, mH = ev.heavy.mass, M = mL + mH;
  double scale = half / std::log10(1 + traj.back().separation);
  for (size_t i = 0; i < traj.size(); ++i) {
    const auto& p = traj[i];
    int cl = int(std::lround(scale * std::log10(1 + p.separation * mH / M)));
    int ch = int(std::lround(scale * std::log10(1 + p.separation * mL / M)));
    std::string line(2 * half + 1, ' ');
    line[half] = '+';
    std::string out;
    for (int c = 0; c <= 2 * half; ++c) {
      if (c == half - std::min(cl, half)) out += "\033[36m●\033[0m";
      else if (c == half + std::min(ch, half)) out += "\033[33m●\033[0m";
      else out += line[c] == '+' ? "┼" : "·";
    }
    std::printf("  %9.3g zs  %s\n", p.time, out.c_str());
  }
  std::printf("             \033[36m●\033[0m %s (light)%*s\033[33m●\033[0m %s (heavy)\n",
              nuclide(ev.light.Z, ev.light.A).c_str(), 34, "",
              nuclide(ev.heavy.Z, ev.heavy.A).c_str());
}

void oneEvent(const Compound& cn, std::mt19937_64& rng, std::vector<TrajectoryPoint>* trajOut) {
  section("2.2  One fission event, step by step");
  FissionEvent ev = simulateFission(cn, rng);
  const Fragment &L = ev.light, &H = ev.heavy;
  std::printf("  %s* → %s + %s\n", nuclide(cn.Z, cn.A).c_str(), nuclide(L.Z, L.A).c_str(),
              nuclide(H.Z, H.A).c_str());
  std::printf("  Q-value (ground states)   %7.2f MeV   mass difference × c²\n", ev.q);
  std::printf("  + excitation E*           %7.2f MeV\n", cn.excitation);
  std::printf("  Scission distance         %7.2f fm    (touching spheres would be %.2f fm)\n",
              ev.scissionDistance, nuclearRadius(L.A) + nuclearRadius(H.A));
  std::printf("  Coulomb energy at scission%7.2f MeV → total kinetic energy (TKE)\n", ev.tke);
  std::printf("  Fragment excitation       %7.2f MeV   (%s %.2f, %s %.2f)\n",
              L.excitation + H.excitation, nuclide(L.Z, L.A).c_str(), L.excitation,
              nuclide(H.Z, H.A).c_str(), H.excitation);

  auto traj = coulombTrajectory(ev, 1e4);
  if (trajOut) *trajOut = traj;
  std::printf("\n  Coulomb acceleration (relativistic two-body, RK4):\n");
  std::printf("       time (zs)     time (s)   separation (fm)   KE light   KE heavy   "
              "converted\n");
  for (size_t i = 0; i < traj.size(); ++i) {
    const auto& p = traj[i];
    std::printf("  %14.4g  %11.3g  %16.4g  %9.2f  %9.2f  %8.2f%%\n", p.time, p.time * 1e-21,
                p.separation, p.kineticLight, p.kineticHeavy,
                100 * (p.kineticLight + p.kineticHeavy) / ev.tke);
  }
  std::printf("  Final speeds: light %.4f c (%.0f km/s), heavy %.4f c (%.0f km/s)\n", L.beta,
              L.beta * 299792.458, H.beta, H.beta * 299792.458);
  filmStrip(ev, traj);

  std::printf("\n  Prompt neutrons (evaporated from the moving fragments, ~10⁻¹⁸–10⁻¹⁴ s):\n");
  for (const Fragment* f : {&L, &H}) {
    for (size_t n = 0; n < f->neutronCmEnergies.size(); ++n)
      std::printf("    from %-6s S_n %5.2f MeV, %5.2f MeV in fragment frame, %5.2f MeV lab\n",
                  nuclide(f->Z, f->A - int(n)).c_str(), f->neutronSeparations[n],
                  f->neutronCmEnergies[n], f->neutronLabEnergies[n]);
    std::printf("    %-6s spins at %dħ (%.2f MeV rotational), sheds %.2f MeV as prompt γ → %s\n",
                nuclide(f->Z, f->A).c_str(), f->spin, f->rotationalEnergy, f->gammaEnergy,
                nuclide(f->Z, f->finalA).c_str());
  }
  std::printf("  %d neutrons released — each can split another nucleus (chain reaction).\n",
              ev.neutrons());

  std::printf("\n  β⁻ decay to stability (seconds to years later; half-lives not modelled):\n");
  for (const Fragment* f : {&L, &H}) {
    std::string chain = "  " + nuclide(f->Z, f->finalA);
    double total = 0;
    for (const auto& b : f->betaChain) {
      char step[64];
      std::snprintf(step, sizeof step, " →(%.1f) %s", b.q, nuclide(b.Z + 1, f->finalA).c_str());
      chain += step;
      total += b.q;
    }
    std::printf("%s   [%.1f MeV]\n", chain.c_str(), total);
  }
}

// Watt spectrum fit to U-235(n_th,f), ENDF: N(E) ∝ exp(-E/a) sinh(sqrt(bE)).
double wattBin(double lo, double hi) {
  auto watt = [](double e) { return std::exp(-e / 0.988) * std::sinh(std::sqrt(2.249 * e)); };
  double sum = 0, norm = 0;
  for (int i = 0; i < 4000; ++i) {
    double e = (i + 0.5) * 20.0 / 4000;
    norm += watt(e);
    if (e >= lo && e < hi) sum += watt(e);
  }
  return sum / norm;
}

RunStats monteCarlo(const Compound& cn, int events, std::uint64_t seed, int threads) {
  if (threads <= 0) threads = int(std::max(1u, std::thread::hardware_concurrency()));
  auto start = std::chrono::steady_clock::now();
  RunStats s = runMonteCarlo(cn, events, seed, threads);
  double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  char title[128];
  std::snprintf(title, sizeof title, "2.3  Monte Carlo: %d fission events (%.2f s, %d threads)",
                events, seconds, threads);
  section(title);
  bool reference = cn.targetZ == 92 && cn.targetA == 235 && cn.neutronEnergy < 1e-3;

  std::printf("  quantity                          model     %s\n",
              reference ? "measured U-235(n_th,f)" : "");
  auto row = [&](const char* name, double model, const char* measured) {
    std::printf("  %-30s %8.2f     %s\n", name, model, reference ? measured : "");
  };
  row("prompt neutrons per fission ν", s.meanNu, "2.41    (ENDF/B-VIII.0)");
  row("total kinetic energy (MeV)", s.meanTke, "170.5   (Madland 2006)");
  row("  light fragment KE (MeV)", s.meanLightKe, "≈100");
  row("  heavy fragment KE (MeV)", s.meanHeavyKe, "≈70");
  row("mean light fragment mass", s.meanLightA, "96.4");
  row("mean heavy fragment mass", s.meanHeavyA, "139.6");
  row("mean neutron lab energy (MeV)", s.meanNeutronLabEnergy, "2.0");
  row("prompt γ energy (MeV)", s.meanPromptGamma, "6.6–6.9");
  row("mean fragment spin (ħ)", s.meanSpin, "≈7     (Wilhelmy 1972)");

  std::printf("\n  Fragment mass yield (pre-neutron, %% per fission, 2-mass bins):\n");
  double peak = 0;
  for (int a = 70; a < 172; a += 2) peak = std::fmax(peak, s.yieldPre[a] + s.yieldPre[a + 1]);
  for (int a = 72; a < 168; a += 2) {
    double y = s.yieldPre[a] + s.yieldPre[a + 1];
    std::printf("  A %3d-%3d │%s %5.2f%%\n", a, a + 1, bar(y / peak, 50).c_str(), 100 * y);
  }

  std::printf("\n  Prompt neutron multiplicity P(ν):\n");
  static const double kMeasuredP[] = {0.032, 0.172, 0.336, 0.303, 0.129, 0.025, 0.003, 0.0};
  for (int nu = 0; nu < 8; ++nu)
    std::printf("  ν=%d │%s %5.3f   %s\n", nu, bar(s.multiplicity[nu] / 0.5, 40).c_str(),
                s.multiplicity[nu],
                reference ? (std::to_string(kMeasuredP[nu]).substr(0, 5) + " measured").c_str()
                          : "");

  std::printf("\n  Prompt neutron lab spectrum, model vs the evaluated Watt spectrum (ENDF):\n");
  for (int b = 0; b < 32; b += 2) {
    double model = s.neutronSpectrum[b] + s.neutronSpectrum[b + 1];
    double watt = wattBin(b * 0.25, (b + 2) * 0.25);
    std::printf("  %4.1f-%-4.1f MeV │%s %5.3f   Watt %5.3f\n", b * 0.25, (b + 2) * 0.25,
                bar(model / 0.25, 40).c_str(), model, watt);
  }

  section("2.4  Where the energy goes (average per fission)");
  double prompt = s.meanTke + s.meanNeutronKe + s.meanPromptGamma;
  std::printf("  Kinetic energy of fragments     %7.2f MeV\n", s.meanTke);
  std::printf("  Prompt neutron kinetic energy   %7.2f MeV\n", s.meanNeutronKe);
  std::printf("  Prompt γ rays                   %7.2f MeV\n", s.meanPromptGamma);
  std::printf("  β⁻ chains (e⁻, ν̄, delayed γ)    %7.2f MeV\n", s.meanBetaChain);
  std::printf("  Total                           %7.2f MeV   (measured ≈ 200 MeV, ~10 of it\n"
              "                                                 carried off by antineutrinos)\n",
              prompt + s.meanBetaChain);
  std::printf("  Check: Q + E* - Σ S_n = %.2f MeV equals the prompt total %.2f MeV\n",
              s.meanQ + cn.excitation - s.meanNeutronBinding, prompt);
  std::printf("  1 g of U-235 fully fissioned ≈ %.2g J ≈ %.0f t of TNT\n",
              (prompt + s.meanBetaChain) * 1.602176634e-13 * 6.02214076e23 / 235.0439,
              (prompt + s.meanBetaChain) * 1.602176634e-13 * 6.02214076e23 / 235.0439 / 4.184e9);
  return s;
}

// Plot data for scripts/plots.py.
void writeCsv(const std::string& dir, const RunStats& s, const std::vector<TrajectoryPoint>& traj) {
  auto open = [&](const char* name) {
    std::string path = dir + "/" + name;
    FILE* f = std::fopen(path.c_str(), "w");
    if (!f) {
      std::fprintf(stderr, "atomsim: cannot write %s\n", path.c_str());
      std::exit(1);
    }
    return f;
  };
  FILE* f = open("yields.csv");
  std::fprintf(f, "A,pre_neutron,post_neutron\n");
  for (size_t a = 60; a < 180; ++a) std::fprintf(f, "%zu,%g,%g\n", a, s.yieldPre[a], s.yieldPost[a]);
  std::fclose(f);

  f = open("multiplicity.csv");
  static const double kMeasuredP[] = {0.032, 0.172, 0.336, 0.303, 0.129, 0.025, 0.003, 0.0};
  std::fprintf(f, "nu,model,measured\n");
  for (int nu = 0; nu < 8; ++nu) std::fprintf(f, "%d,%g,%g\n", nu, s.multiplicity[nu], kMeasuredP[nu]);
  std::fclose(f);

  f = open("spectrum.csv");
  std::fprintf(f, "E_low,E_high,model,watt\n");
  for (size_t b = 0; b < 40; ++b)
    std::fprintf(f, "%g,%g,%g,%g\n", b * 0.25, (b + 1) * 0.25, s.neutronSpectrum[b],
                 wattBin(b * 0.25, (b + 1) * 0.25));
  std::fclose(f);

  f = open("trajectory.csv");
  std::fprintf(f, "time_zs,separation_fm,ke_light,ke_heavy,potential\n");
  for (const auto& p : traj)
    std::fprintf(f, "%g,%g,%g,%g,%g\n", p.time, p.separation, p.kineticLight, p.kineticHeavy,
                 p.potential);
  std::fclose(f);

  f = open("hydrogen.csv");
  std::vector<RadialState> states;
  for (auto [n, l] : {std::pair{1, 0}, {2, 0}, {2, 1}, {3, 0}, {3, 2}})
    states.push_back(solveHydrogenic(n, l, 1.0, phys::kHydrogenReducedMass));
  std::fprintf(f, "r_bohr,1s,2s,2p,3s,3d\n");
  for (double r = 0; r <= 30; r += 0.05) {
    std::fprintf(f, "%g", r);
    for (const auto& st : states) {
      size_t i = size_t(r / st.step);
      std::fprintf(f, ",%g", i < st.u.size() ? st.u[i] * st.u[i] : 0.0);
    }
    std::fprintf(f, "\n");
  }
  std::fclose(f);
}

}  // namespace

int main(int argc, char** argv) {
  useDataDirectoryIfPresent((std::filesystem::path(argv[0]).parent_path() / "data").string());
  Options o = parse(argc, argv);
  std::mt19937_64 rng(o.seed);
  Compound cn = formCompound(o.targetZ, o.targetA, o.neutronEnergy);

  heading("PART 1 · THE ATOM");
  hydrogenAtom();
  heavyAtom(o.targetZ, o.targetA);

  heading("PART 2 · SPLITTING THE NUCLEUS");
  capture(cn);
  if (cn.transmission < 1e-3) {
    std::printf(
        "\n  E* is %.2f MeV below the barrier: the nucleus does not split. It de-excites by\n"
        "  γ emission (radiative capture). Try --energy 2 for fast neutrons.\n",
        cn.barrier - cn.excitation);
    return 0;
  }
  std::vector<TrajectoryPoint> traj;
  oneEvent(cn, rng, &traj);
  RunStats stats = monteCarlo(cn, o.events, o.seed, o.threads);
  if (!o.csvDir.empty()) writeCsv(o.csvDir, stats, traj);

  heading("MODEL LIMITS");
  std::puts(
      "  Exact: hydrogen bound states, Coulomb acceleration, conservation laws.\n"
      "  Measured: all nuclear masses (AME2020). Semi-empirical: fission barrier,\n"
      "  mass-yield modes and scission distance (fitted to data). Not modelled: QED,\n"
      "  pre-scission dynamics (needs time-dependent DFT), delayed neutrons, half-lives.");
  return 0;
}
