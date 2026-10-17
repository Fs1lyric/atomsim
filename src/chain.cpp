#include "chain.hpp"

#include "constants.hpp"

#include <algorithm>
#include <cmath>

namespace {

constexpr double kAvogadro = 6.02214076e23;
constexpr double kBarn = 1e-24;  // cm^2

void isotropic(std::mt19937_64& rng, double& dx, double& dy, double& dz) {
  std::uniform_real_distribution<double> u(0.0, 1.0);
  double cosT = 2 * u(rng) - 1, phi = 2 * phys::kPi * u(rng);
  double sinT = std::sqrt(1 - cosT * cosT);
  dx = sinT * std::cos(phi);
  dy = sinT * std::sin(phi);
  dz = cosT;
}

// Integer neutron count with the right mean: floor(nu) or floor(nu) + 1.
int sampleNu(double nu, std::mt19937_64& rng) {
  int base = int(nu);
  return base + (std::uniform_real_distribution<double>(0, 1)(rng) < nu - base ? 1 : 0);
}

double sampleFlight(double sigmaTotal, std::mt19937_64& rng) {
  return -std::log(1.0 - std::uniform_real_distribution<double>(0, 1)(rng)) / sigmaTotal;
}

// Distance from (x,y,z) along unit (dx,dy,dz) to a sphere of radius R around the origin.
double exitDistance(double x, double y, double z, double dx, double dy, double dz, double R) {
  double b = x * dx + y * dy + z * dz, c = x * x + y * y + z * z - R * R;
  return -b + std::sqrt(std::max(0.0, b * b - c));
}

}  // namespace

Material uranium235() {
  // 2 MeV neutrons: v = c sqrt(1 - (m / (m + T))^2).
  double gamma = 1 + 2.0 / phys::kMassNeutron;
  double speed = 2.99792458e10 * std::sqrt(1 - 1 / (gamma * gamma));
  return {18.71, 235.0439, 1.235, 4.566, 0.0, 2.637, speed};
}

double numberDensity(const Material& m) { return m.density * kAvogadro / m.molarMass; }

double sphereRadius(const Material& m, double massKg) {
  return std::cbrt(3 * massKg * 1000 / (4 * phys::kPi * m.density));
}

double multiplicationFactor(const Material& m, double massKg, int perGeneration,
                            std::uint64_t seed) {
  const double R = sphereRadius(m, massKg);
  const double n = numberDensity(m) * kBarn;
  const double sf = n * m.sigmaFission, ss = n * m.sigmaScatter, sc = n * m.sigmaCapture;
  const double st = sf + ss + sc;
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> u(0.0, 1.0);

  // Start uniformly in the sphere; the fission source converges within a few generations.
  struct Site { double x, y, z; };
  std::vector<Site> source(perGeneration);
  for (auto& s : source) {
    double r = R * std::cbrt(u(rng)), dx, dy, dz;
    isotropic(rng, dx, dy, dz);
    s = {r * dx, r * dy, r * dz};
  }

  constexpr int kInactive = 10, kActive = 40;
  double kSum = 0;
  std::vector<Site> bank;
  for (int gen = 0; gen < kInactive + kActive; ++gen) {
    bank.clear();
    double produced = 0;
    for (const Site& s : source) {
      double x = s.x, y = s.y, z = s.z, dx, dy, dz;
      isotropic(rng, dx, dy, dz);
      for (;;) {
        double d = sampleFlight(st, rng);
        if (d >= exitDistance(x, y, z, dx, dy, dz, R)) break;  // leaks out
        x += d * dx;
        y += d * dy;
        z += d * dz;
        double pick = u(rng) * st;
        if (pick < sf) {
          produced += m.nu;
          for (int k = sampleNu(m.nu, rng); k > 0; --k) bank.push_back({x, y, z});
          break;
        }
        if (pick < sf + sc) break;  // captured
        isotropic(rng, dx, dy, dz);  // scattered
      }
    }
    if (gen >= kInactive) kSum += produced / source.size();
    if (bank.empty()) return gen >= kInactive ? kSum / (gen - kInactive + 1) : 0.0;
    // Next generation: perGeneration sites drawn from the fission bank.
    for (auto& s : source) s = bank[std::size_t(u(rng) * bank.size()) % bank.size()];
  }
  return kSum / kActive;
}

double criticalMass(const Material& m, std::uint64_t seed) {
  double lo = 5, hi = 200;
  for (int i = 0; i < 20; ++i) {
    double mid = std::sqrt(lo * hi);
    (multiplicationFactor(m, mid, 20000, seed + i) < 1 ? lo : hi) = mid;
  }
  return std::sqrt(lo * hi);
}

ChainReaction::ChainReaction(const Material& m, double massKg, std::uint64_t seed, int cap)
    : material_(m), radius_(sphereRadius(m, massKg)), cap_(cap), rng_(seed) {
  double n = numberDensity(m) * kBarn;
  sigmaTotal_ = n * (m.sigmaFission + m.sigmaScatter + m.sigmaCapture);
}

Neutron ChainReaction::spawn(double x, double y, double z) {
  Neutron n{x, y, z, 0, 0, 0, sampleFlight(sigmaTotal_, rng_)};
  isotropic(rng_, n.dx, n.dy, n.dz);
  return n;
}

double ChainReaction::distanceToSurface(const Neutron& n) const {
  return exitDistance(n.x, n.y, n.z, n.dx, n.dy, n.dz, radius_);
}

void ChainReaction::inject(int count) {
  for (int i = 0; i < count; ++i) neutrons_.push_back(spawn(0, 0, 0));
}

void ChainReaction::advance(double seconds) {
  std::uniform_real_distribution<double> u(0.0, 1.0);
  const double sf = material_.sigmaFission, sc = material_.sigmaCapture;
  const double st = sf + material_.sigmaScatter + sc;
  flashes_.clear();
  std::vector<Neutron> next;
  next.reserve(neutrons_.size() * 2);
  for (Neutron n : neutrons_) {
    double travel = material_.neutronSpeed * seconds;
    bool alive = true;
    while (alive) {
      double step = std::min(travel, n.toCollision);
      if (distanceToSurface(n) <= step) {  // escapes the sphere
        alive = false;
        break;
      }
      n.x += step * n.dx;
      n.y += step * n.dy;
      n.z += step * n.dz;
      travel -= step;
      n.toCollision -= step;
      if (n.toCollision > 0) break;  // still in flight at the end of the step
      double pick = u(rng_) * st;
      if (pick < sf) {
        flashes_.push_back({n.x, n.y, n.z});
        fissions_ += weight_;
        for (int k = sampleNu(material_.nu, rng_); k > 0; --k) next.push_back(spawn(n.x, n.y, n.z));
        alive = false;
      } else if (pick < sf + sc) {
        alive = false;
      } else {
        isotropic(rng_, n.dx, n.dy, n.dz);
        n.toCollision = sampleFlight(sigmaTotal_, rng_);
      }
    }
    if (alive) next.push_back(n);
  }
  neutrons_.swap(next);
  while (int(neutrons_.size()) > cap_) {  // keep the work bounded; weights keep counts right
    std::shuffle(neutrons_.begin(), neutrons_.end(), rng_);
    neutrons_.resize(neutrons_.size() / 2);
    weight_ *= 2;
  }
  time_ += seconds;
}
