#include "three_fluid.hpp"

void tidal_bench_AB(void){

  {
    ThreeFluidSim sim;
    sim.initSolver(150);
    sim.initCoeffsYiming();

    // Manually set DM mass
    sim.md = 0.1 * sim.ms;

    sim.initPlummerYiming(0.5, 1e-10, 1.0, 1.0, 1.0);

    LagrangianRadiiObserver observer({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    
    std::string dir = "output/baseline_AB/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(50000000, observer);
    observer.save(dir);
  }

  {
    ThreeFluidSim sim;
    sim.initSolver(150);
    sim.initCoeffsYiming();

    // Manually set DM mass
    sim.md = 0.1 * sim.ms;
    sim.tidal_cutoff = 1;
    sim.tidal_cutoff_factor = 1e-3;
    sim.tidal_radius = 5e1;

    sim.initPlummerYiming(0.5, 1e-10, 1.0, 1.0, 1.0);

    LagrangianRadiiObserver observer({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    
    std::string dir = "output/with_tidal_AB/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(50000000, observer);
    observer.save(dir);
  }

}

int main() {
  tidal_bench_AB();
  return 0;
  
  // sim.printParams();
  // sim.initCoeffs();
  
  //ApproximateCentralDensityObserver observer({0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12});
  
  // const std::vector times_to_save({0.0, 5.5, 5.565, 5.56517});
  const std::vector times_to_save({0.0, 5.5, 5.53, 5.538});

  {
    ThreeFluidSim sim;
    sim.initSolver(150);
    sim.initCoeffsYiming();

    // Manually make DM the same as single star
    sim.md = sim.ms;
    sim.c2[FD] = sim.c2[FS];

    sim.initPlummerYiming(1.0, 1e-10, 1e-10, 1.0, 1.0);
    //sim.printParams();

    ApproximateTimeObserver observer1(times_to_save);
    LagrangianRadiiObserver observer2({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    ObserverPack observer(observer1, observer2);
    
    std::string dir = "output/one_fluid_Yiming/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(50000000, observer);
    observer.save(dir);
  }

  // {
  //   ThreeFluidSim sim;
  //   sim.initSolver(150);
  //   sim.initCoeffsYiming();

  //   // Manually make DM the same as single star
  //   sim.md = sim.ms;
  //   sim.c2[FD] = sim.c2[FS];
  
  //   sim.initPlummerYiming(0.5, 1e-10, 1.0, 1.0, 1.0);
    
  //   std::string dir = "output/one_fluid_split_in_two_Yiming/";
  //   prepare_directory_for_output(dir);
  //   sim.saveParams(dir);
  //   ApproximateTimeObserver observer(times_to_save);
  //   sim.evolve(50000000, observer);
  //   observer.save(dir);
  // }


  {
    ThreeFluidSim sim;
    sim.initSolver(150);
    sim.initCoeffsYiming();

    // Manually make DM the same as single star
    sim.md = sim.ms;
    sim.c2[FD] = sim.c2[FS];
    sim.tidal_cutoff = 1;

    // sim.printParams();
    sim.initPlummerYiming(1.0, 1e-10, 1e-10, 1.0, 1.0);

    ApproximateTimeObserver observer1(times_to_save);
    LagrangianRadiiObserver observer2({0.01, 0.05, 0.1, 0.2, 0.5, 0.7});
    ObserverPack observer(observer1, observer2);
    
    std::string dir = "output/tidal_effect_Yiming/";
    prepare_directory_for_output(dir);
    sim.saveParams(dir);
    sim.evolve(100000, observer);
    observer.save(dir);
  }


  return 0;
}
