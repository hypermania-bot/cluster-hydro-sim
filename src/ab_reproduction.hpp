#pragma once
#include "three_fluid.hpp"
#include <cmath>
#include <numbers>
#include <stdexcept>

// Physical assumptions of the AB comparison, not dimensionful solver state.
struct ABInitParam {
  double particle_number=50000; // Reference untruncated cluster; finite domain retains less.
  double stellar_mass_msun=1, dm_mass_msun=0.1;
  double half_mass_radius_pc=3;
  double galactic_radius_pc=8500, circular_speed_kms=220;
  double outer_radius_over_initial_jacobi=0.45;
};

struct ABObserverParam {
  double length_unit_pc=1, time_unit_myr=1, friction_time_myr=1;
  double radius_to_nbody=3*std::numbers::pi/16;
  std::array<double,6> mass_fractions={0.01,0.05,0.1,0.2,0.5,0.7};
};

inline ABObserverParam initializeAB(ThreeFluidSim& sim,const ABInitParam& init,bool tide) {
  for(double x:{init.particle_number,init.stellar_mass_msun,init.dm_mass_msun,
                init.half_mass_radius_pc,init.galactic_radius_pc,init.circular_speed_kms,
                init.outer_radius_over_initial_jacobi})
    if(!std::isfinite(x)||x<=0)throw std::invalid_argument("invalid AB initialization parameter");
  if(init.outer_radius_over_initial_jacobi>=1)
    throw std::invalid_argument("AB boundary must lie inside the initial Jacobi radius");
  constexpr double G=4.498502151575286e-3; // pc^3 / Msun / Myr^2
  constexpr double kms_to_pc_myr=1.022712165045695;
  const double mass=2*init.particle_number/
    (1/init.stellar_mass_msun+1/init.dm_mass_msun); // equal component masses
  const double mass_ratio=init.dm_mass_msun/init.stellar_mass_msun;
  ABObserverParam observing;
  observing.length_unit_pc=init.half_mass_radius_pc*std::sqrt(std::pow(2.,2./3)-1);
  const double a=observing.length_unit_pc;
  const double lnL=std::log(0.8*mass/init.stellar_mass_msun/(1+mass_ratio));
  if(lnL<=0)throw std::invalid_argument("nonpositive AB Coulomb logarithm");
  observing.time_unit_myr=std::sqrt(mass/(3*init.stellar_mass_msun))*std::pow(a,1.5)/
    (std::sqrt(G*init.stellar_mass_msun)*lnL);
  observing.friction_time_myr=0.035*std::pow(init.half_mass_radius_pc,1.5)*
    std::sqrt(mass/G)/init.stellar_mass_msun; // AB Eq. (2)
  const double physical_q=2*std::pow(init.circular_speed_kms*kms_to_pc_myr/init.galactic_radius_pc,2);
  const double q=physical_q*std::pow(a,3)/(3*G*mass);
  const double jacobi=std::sqrt(std::pow(1/(3*q),2./3)-1); // continuum Plummer
  // Use the same finite domain in the q=0 control; not a mass-removal law.
  const double outer=init.outer_radius_over_initial_jacobi*jacobi;
  sim.param.tidal_q=tide?q:0;
  sim.initSolver(sim.param.N);
  sim.initPlummer(0.5,1e-10,1,1,1,outer);
  // Fixed continuum mass unit M0=3*M, not the evolving/truncated mass.
  sim.param.ms=init.stellar_mass_msun/(3*mass);
  sim.param.mb=2*sim.param.ms;sim.param.md=mass_ratio*sim.param.ms;
  sim.initCoeffs(mass/init.stellar_mass_msun,2,mass_ratio);
  sim.param.c4.fill(0);
  return observing;
}
