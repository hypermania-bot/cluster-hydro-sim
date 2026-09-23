#pragma once
#include "moving_initialization.hpp"
#include <filesystem>
#include <fstream>

struct MovingComparisonInit {
  long long sample=0; // 0 single Yiming, 1 split, 2 AB no stripping, 3 canonical
  long long canonical_zones=500;
  double stars=1e6;
};
struct MovingComparisonObserverParam {double interval=0.02;long long steps_between_snapshots=1000;};

inline void initializeMovingComparison(ThreeFluidSim& s,const MovingComparisonInit& p) {
  if(p.sample<0||p.sample>3)throw std::invalid_argument("comparison sample must be 0..3");
  s.initSolver(p.sample==3?p.canonical_zones:150);
  if(p.sample<3) {
    s.initCoeffsYiming();s.param.md=(p.sample==2?0.1:1)*s.param.ms;
    if(p.sample<2)s.param.c2[FD]=s.param.c2[FS];
    s.initPlummerYiming(p.sample==0?1:0.5,1e-10,p.sample==0?1e-10:1,1,1);
  } else {
    s.initPlummer(1,1e-10,1e-10,1,1);
    double total=0;for(int f=0;f<NF;++f)total+=s.Menc[f][s.param.N-1];
    s.param.ms=total/p.stars;s.param.mb=2*s.param.ms;s.param.md=1e-10*s.param.ms;
    s.initCoeffs(p.stars,2,1e-10);
  }
}

// Both overloads write the same profile schema; radii are CELL CENTRES,
// sigma is later plotted as sqrt(2U/3), and times are actual accepted times.
struct MovingComparisonObserver {
  MovingComparisonObserverParam param;
  std::string directory;
  std::ofstream manifest, history;
  double next_time=0;
  long long count=0;
  std::vector<double> latest;
  double latest_time=0,last_saved_time=-1;
  int latest_zones=0;
  MovingComparisonObserver(const std::string& dir,MovingComparisonObserverParam p):param(p),directory(dir) {
    std::filesystem::create_directories(dir);
    manifest.open(dir+"/snapshots.csv");history.open(dir+"/history.csv");
    manifest.exceptions(std::ios::badbit|std::ios::failbit);
    history.exceptions(std::ios::badbit|std::ios::failbit);
    manifest<<std::setprecision(17)<<"snapshot,time,zones\n";
    history<<std::setprecision(17)<<"step,time,dt,rho_s,u_s,rho_b,u_b,rho_d,u_d,m_s,m_b,m_d,out_s,out_b,out_d,linear_residual,last_change,limiting_index,energy_ledger_error,gravity_work,heating\n";
  }
  void record(double time,int n,const std::vector<double>& data) {
    if(time+1e-14<next_time)return;
    std::ofstream out(directory+"/profile_"+std::to_string(count)+".dat",std::ios::binary);
    out.exceptions(std::ios::badbit|std::ios::failbit);
    out.write(reinterpret_cast<const char*>(data.data()),data.size()*sizeof(double));
    manifest<<count++<<','<<time<<','<<n<<'\n';manifest.flush();
    last_saved_time=time;
    next_time=time+param.interval;
  }
  void saveFinal() {
    if(latest_zones>0&&latest_time>last_saved_time) {next_time=0;record(latest_time,latest_zones,latest);}
    history.flush();
  }
  void operator()(const ThreeFluidSim& s) {
    history<<s.step<<','<<s.totalTime<<','<<s.Deltat;
    for(int f=0;f<NF;++f)history<<','<<s.Rho[f][0]<<','<<s.U[f][0];
    for(int f=0;f<NF;++f)history<<','<<s.Menc[f][s.param.N-1];
    history<<",0,0,0,0,0,0,0,0,0\n";
    latest.clear();latest.reserve(15*s.param.N);
    for(int i=0;i<s.param.N;++i)for(int f=0;f<NF;++f) {
      latest.insert(latest.end(),{(s.R[f][i]+(i?s.R[f][i-1]:0))/2,s.Rho[f][i],s.U[f][i],s.Menc[f][i],0});
    }
    latest_time=s.totalTime;latest_zones=s.param.N;record(latest_time,latest_zones,latest);
  }
  void operator()(const MovingThreeFluidSim& s) {
    history<<s.step<<','<<s.totalTime<<','<<s.Deltat;
    for(int f=0;f<NF;++f)history<<','<<s.value(0,f,s.RHO)<<','<<s.value(0,f,s.U);
    for(int f=0;f<NF;++f)history<<','<<s.value(s.zones()-1,f,s.MASS);
    for(int f=0;f<NF;++f)history<<','<<s.escaped_mass[f];
    history<<','<<s.linear_residual<<','<<s.last_change<<','<<s.limiting_index
           <<','<<s.energy_ledger_error<<','<<s.last_gravity_work<<','<<s.last_heating<<'\n';
    latest.clear();latest.reserve(15*s.zones());
    for(int i=0;i<s.zones();++i)for(int f=0;f<NF;++f)
      latest.insert(latest.end(),{s.radius(i),s.value(i,f,s.RHO),s.value(i,f,s.U),s.value(i,f,s.MASS),s.value(i,f,s.VEL)});
    latest_time=s.totalTime;latest_zones=s.zones();
    if(s.step%param.steps_between_snapshots==0)next_time=0;
    record(latest_time,latest_zones,latest);
  }
};

inline void runMovingComparison(const std::string& directory,int sample,double until,
                               double epsilon_override=0,double max_dt=1e-3,int zones=500,
                               long long max_steps=20000) {
  MovingComparisonInit initial;initial.sample=sample;initial.canonical_zones=zones;
  if(!(std::isfinite(until)&&until>0&&std::isfinite(max_dt)&&max_dt>0&&
       std::isfinite(epsilon_override)&&epsilon_override>=0&&zones>=3&&max_steps>0))
    throw std::invalid_argument("invalid moving comparison controls");
  if(std::filesystem::exists(directory+"/hydro/snapshots.csv")||
     std::filesystem::exists(directory+"/moving/snapshots.csv"))
    throw std::invalid_argument("comparison output already exists; choose a fresh directory");
  ThreeFluidSim hydro;initializeMovingComparison(hydro,initial);
  hydro.param.maxTime=until;hydro.param.max_timestep=max_dt;
  hydro.param.Deltat=std::min(1e-5,max_dt);hydro.param.u_change_tolerance=1e-3;
  MovingThreeFluidSim moving;
  // t0 = r0/sqrt(u0)/(3 m_s_hat ln Lambda_sd), with m_s in M0 units.
  // Use the actual stored mass/count, also for the legacy Yiming preset;
  // initial.stars sets the count only in the canonical initializer.
  double total_mass=0;for(int f=0;f<NF;++f)total_mass+=hydro.Menc[f][hydro.param.N-1];
  const double lnLambda=std::log(GAMMA_LAMBDA*2*total_mass/(hydro.param.ms+hydro.param.md));
  if(!(lnLambda>0))throw std::invalid_argument("nonpositive reference Coulomb logarithm");
  moving.param.epsilon=epsilon_override>0?epsilon_override:std::pow(3*hydro.param.ms*lnLambda,2);
  moving.param.Deltat=hydro.param.Deltat;moving.param.max_timestep=max_dt;
  moving.param.maxSteps=max_steps;hydro.param.maxSteps=max_steps;
  moving.param.maxTime=until;initializeMovingFromHydrostatic(moving,hydro);
  const MovingComparisonObserverParam observing;
  std::filesystem::create_directories(directory+"/initialization");
  std::filesystem::create_directories(directory+"/observer");
  std::filesystem::create_directories(directory+"/hydro");
  std::filesystem::create_directories(directory+"/moving");
  save_param_for_Mathematica(initial,directory+"/initialization/");
  save_param_for_Mathematica(observing,directory+"/observer/");
  hydro.saveParams(directory+"/hydro/");
  save_param_for_Mathematica(moving.param,directory+"/moving/");
  {
    MovingComparisonObserver observer(directory+"/hydro",observing);
    hydro.evolve(observer);
    observer.saveFinal();
    std::cout<<"hydro sample="<<sample<<" time="<<hydro.totalTime<<" steps="<<hydro.step<<std::endl;
    if(hydro.totalTime<until)throw std::runtime_error("hydro comparison did not reach requested time");
  }
  {
    MovingComparisonObserver observer(directory+"/moving",observing);
    try {moving.evolve(observer);}
    catch(const std::exception& e) {
      observer.saveFinal();
      std::cerr<<"moving failed at t="<<std::setprecision(17)<<moving.totalTime
               <<" dt="<<moving.Deltat<<" step="<<moving.step<<": "<<e.what()<<'\n';throw;
    }
    observer.saveFinal();
    std::cout<<"moving sample="<<sample<<" time="<<moving.totalTime<<" steps="<<moving.step<<std::endl;
    if(moving.totalTime<until)throw std::runtime_error("moving comparison did not reach requested time");
  }
}
