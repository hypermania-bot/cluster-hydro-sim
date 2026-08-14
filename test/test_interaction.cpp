#include "../src/three_fluid.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

int main() {
  ThreeFluidSim sim;
  sim.initSolver(5);
  sim.Deltat = 1e-3;
  sim.ms = 1.0;
  sim.md = 0.1;

  for (int f = 0; f < NF; ++f) {
    sim.c2[f] = 0.0;
    for (int f2 = 0; f2 < NF; ++f2) {
      sim.c1[f][f2] = 0.0;
      sim.c4[f][f2] = 0.0;
    }
  }
  sim.c1[FS][FD] = sim.c1[FD][FS] = C1;

  for (int f = 0; f < NF; ++f) {
    sim.R[f] = Eigen::VectorXd::LinSpaced(sim.N, 1.0, 5.0);
    sim.U[f] = Eigen::VectorXd::Ones(sim.N);
    sim.Rho[f] = Eigen::VectorXd::Ones(sim.N);
  }
  sim.Rho[FS].setConstant(2.0);
  sim.Rho[FD].setConstant(3.0);

  constexpr int zone = 2;
  const double energyBefore =
      sim.Rho[FS][zone] * sim.U[FS][zone]
      + sim.Rho[FD][zone] * sim.U[FD][zone];

  sim.solveConductionLAPACKE();

  const double energyAfter =
      sim.Rho[FS][zone] * sim.U[FS][zone]
      + sim.Rho[FD][zone] * sim.U[FD][zone];
  const double relativeError = std::abs(energyAfter - energyBefore) / energyBefore;

  if (!(sim.U[FS][zone] < 1.0 && sim.U[FD][zone] > 1.0)) {
    std::cerr << "Interaction must cool the heavy fluid and heat the light fluid\n";
    return 1;
  }
  if (relativeError > 1e-12) {
    std::cerr << "Interaction energy conservation error: " << relativeError << '\n';
    return 1;
  }

  std::cout << "interaction regression passed; relative energy error = "
            << relativeError << '\n';

  ThreeFluidSim remap;
  remap.initSolver(20);
  remap.initPlummer();
  const double massBefore[NF] = {
      remap.Menc[FS][remap.N - 1],
      remap.Menc[FB][remap.N - 1],
      remap.Menc[FD][remap.N - 1]};
  remap.R[FS] *= 0.98;
  remap.R[FD] *= 1.05;
  remap.realign();

  for (int f = 0; f < NF; ++f) {
    const double massError =
        std::abs(remap.Menc[f][remap.N - 1] - massBefore[f]) / massBefore[f];
    if (massError > 1e-12) {
      std::cerr << "Realignment mass conservation error for fluid " << f
                << ": " << massError << '\n';
      return 1;
    }
  }
  std::cout << "realignment regression passed\n";
  return 0;
}
