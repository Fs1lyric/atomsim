#pragma once

// CODATA 2018 recommended values. Energies in MeV, lengths in fm, times in zs
// (1 zeptosecond = 1e-21 s) unless a name says otherwise.
namespace phys {

inline constexpr double kPi = 3.14159265358979323846;
inline constexpr double kHbarC = 197.3269804;               // MeV fm
inline constexpr double kAlpha = 7.2973525693e-3;           // fine-structure constant
inline constexpr double kCoulomb = kAlpha * kHbarC;         // e^2/(4 pi eps0) = 1.43996 MeV fm
inline constexpr double kSpeedOfLight = 299.792458;         // fm per zs
inline constexpr double kPlanckEVs = 4.135667696e-15;       // eV s
inline constexpr double kHartreeEV = 27.211386245988;       // eV
inline constexpr double kBohrRadiusFm = 52917.721090;       // fm

inline constexpr double kAtomicMassUnit = 931.49410242;     // MeV/c^2
inline constexpr double kMassNeutron = 939.56542052;        // MeV/c^2
inline constexpr double kMassProton = 938.27208816;         // MeV/c^2
inline constexpr double kMassElectron = 0.51099895000;      // MeV/c^2

// Mass excesses (AME2020) of the free neutron and the hydrogen atom.
inline constexpr double kExcessNeutron = 8.0713181;         // MeV
inline constexpr double kExcessHydrogen = 7.2889710;        // MeV
inline constexpr double kNeutronMinusHydrogen = kExcessNeutron - kExcessHydrogen;  // 0.78235 MeV

// Electron-proton reduced mass, in MeV and in electron masses.
inline constexpr double kHydrogenReducedMassMeV =
    kMassElectron * kMassProton / (kMassElectron + kMassProton);
inline constexpr double kHydrogenReducedMass = kHydrogenReducedMassMeV / kMassElectron;

}  // namespace phys
