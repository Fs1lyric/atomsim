#include "hydrogen.hpp"

#include "constants.hpp"

#include <cmath>
#include <stdexcept>

namespace {

struct Problem {
  int l;
  double Z, mu, E, h;
  int steps;

  double g(double r) const { return l * (l + 1) / (r * r) - 2 * mu * Z / r - 2 * mu * E; }
};

// Numerov integration outward from r = 0. Returns the number of sign changes
// of u on (0, r_max); stores u when out is given.
int integrate(const Problem& p, std::vector<double>* out) {
  const double h = p.h, k = h * h / 12;
  // Numerov works on y = (1 - h^2 g / 12) u. At r = 0, u ~ r^(l+1), so u = 0
  // but g*u has a finite limit: -2 mu Z for l = 0, l(l+1) = 2 for l = 1.
  double limGu = p.l == 0 ? -2 * p.mu * p.Z : (p.l == 1 ? 2.0 : 0.0);
  double yPrev = -k * limGu;
  double u = std::pow(h, p.l + 1) * (1 - p.mu * p.Z * h / (p.l + 1));
  double y = (1 - k * p.g(h)) * u;

  if (out) {
    out->assign(p.steps + 1, 0.0);
    (*out)[1] = u;
  }
  int nodes = 0;
  for (int i = 1; i < p.steps; ++i) {
    double r = i * h;
    double yNext = 2 * y - yPrev + h * h * p.g(r) * u;
    double uNext = yNext / (1 - k * p.g(r + h));
    if (uNext * u < 0) ++nodes;
    yPrev = y;
    y = yNext;
    u = uNext;
    if (out) {
      (*out)[i + 1] = u;
    } else if (std::fabs(u) > 1e200) {  // only the sign matters while shooting
      yPrev *= 1e-200;
      y *= 1e-200;
      u *= 1e-200;
    }
  }
  return nodes;
}

}  // namespace

RadialState solveHydrogenic(int n, int l, double Z, double reducedMass) {
  if (l < 0 || l > 2 || l >= n) throw std::invalid_argument("need 0 <= l <= min(2, n-1)");
  Problem p{l, Z, reducedMass, 0.0, 1e-3, 0};
  // Far enough that e^{-r/n} has died away, close enough that round-off in E
  // does not blow up the growing solution.
  double rMax = 2.0 * n * (n + 8) / (reducedMass * Z);
  p.steps = static_cast<int>(rMax / p.h);

  // A bound state with n - l - 1 radial nodes. Below the eigenvalue u has that
  // many nodes; just above it, one more appears near r_max. Bisect on that.
  const int radialNodes = n - l - 1;
  double lo = -reducedMass * Z * Z, hi = 0.0;
  for (int it = 0; it < 200 && hi - lo > 1e-15 * std::fabs(lo); ++it) {
    p.E = 0.5 * (lo + hi);
    (integrate(p, nullptr) > radialNodes ? hi : lo) = p.E;
  }
  p.E = 0.5 * (lo + hi);

  RadialState s;
  s.n = n;
  s.l = l;
  s.step = p.h;
  s.energyEV = p.E * phys::kHartreeEV;
  integrate(p, &s.u);

  double norm = 0, r1 = 0;
  for (size_t i = 0; i < s.u.size(); ++i) {
    double w = (i == 0 || i + 1 == s.u.size()) ? 0.5 : 1.0;  // trapezoid
    norm += w * s.u[i] * s.u[i];
    r1 += w * i * p.h * s.u[i] * s.u[i];
  }
  norm *= p.h;
  r1 *= p.h;
  for (double& v : s.u) v /= std::sqrt(norm);
  s.meanRadius = r1 / norm;
  return s;
}

double schrodingerEnergyEV(int n, double Z, double reducedMass) {
  return -reducedMass * Z * Z * phys::kHartreeEV / (2.0 * n * n);
}

double diracEnergyEV(int n, double j, double Z, double reducedMassMeV) {
  double k = j + 0.5, za = Z * phys::kAlpha;
  double d = n - k + std::sqrt(k * k - za * za);
  return reducedMassMeV * 1e6 * (1.0 / std::sqrt(1.0 + za * za / (d * d)) - 1.0);
}
