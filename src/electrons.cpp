#include "electrons.hpp"

#include <algorithm>
#include <vector>

namespace {

struct Subshell {
  int n, l, electrons;
};

constexpr char kLetters[] = "spdf";

// Madelung order: increasing n + l, ties broken by increasing n.
std::vector<Subshell> madelungOrder() {
  std::vector<Subshell> order;
  for (int sum = 1; sum <= 8; ++sum)
    for (int n = 1; n <= sum; ++n) {
      int l = sum - n;
      if (l < n && l <= 3) order.push_back({n, l, 0});
    }
  std::stable_sort(order.begin(), order.end(), [](const Subshell& a, const Subshell& b) {
    return a.n + a.l < b.n + b.l || (a.n + a.l == b.n + b.l && a.n < b.n);
  });
  return order;
}

}  // namespace

std::string electronConfiguration(int Z) {
  std::vector<Subshell> shells = madelungOrder();
  int left = Z;
  for (auto& s : shells) {
    s.electrons = std::min(left, 2 * (2 * s.l + 1));
    left -= s.electrons;
  }

  // Measured exceptions: in U and Np a 5f electron sits in 6d instead.
  if (Z == 92 || Z == 93)
    for (auto& s : shells) {
      if (s.n == 5 && s.l == 3) --s.electrons;
      if (s.n == 6 && s.l == 2) ++s.electrons;
    }

  static const int kNobleZ[] = {86, 54, 36, 18, 10, 2};
  static const char* kNobleName[] = {"Rn", "Xe", "Kr", "Ar", "Ne", "He"};
  std::string out;
  size_t first = 0;
  for (int i = 0; i < 6; ++i) {
    if (kNobleZ[i] >= Z) continue;
    int count = 0;
    while (count < kNobleZ[i]) count += 2 * (2 * shells[first++].l + 1);
    out = std::string("[") + kNobleName[i] + "]";
    break;
  }

  std::vector<Subshell> valence(shells.begin() + first, shells.end());
  std::sort(valence.begin(), valence.end(), [](const Subshell& a, const Subshell& b) {
    return a.n < b.n || (a.n == b.n && a.l < b.l);
  });
  for (const auto& s : valence) {
    if (s.electrons == 0) continue;
    if (!out.empty()) out += ' ';
    out += std::to_string(s.n) + kLetters[s.l] + std::to_string(s.electrons);
  }
  return out;
}
