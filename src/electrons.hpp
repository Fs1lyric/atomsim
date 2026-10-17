#pragma once

#include <string>

// Ground-state electron configuration, e.g. "[Rn] 5f3 6d1 7s2".
// Filled by the Madelung (n + l) rule, with the measured exceptions for the
// actinides this program uses.
std::string electronConfiguration(int Z);
