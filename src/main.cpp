#include "three_fluid.hpp"

int main() {
    
  ThreeFluidSim sim("output/split_in_two_Yiming/");
  sim.initSolver(150);
  sim.initCoeffsYiming();
  sim.initPlummerYiming(1e-6, 1e-10, 1.0, 1.0);
  // sim.initPlummerYiming(1e-1, 1.0, 1.0, 1.0);
  sim.updateEnclosedMass();
  sim.printParams();
  // sim.initCoeffs();

  sim.saveParams();

  ApproximateCentralDensityObserver observer({0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12});
    
  sim.evolve(50000000, observer);

  return 0;
}
