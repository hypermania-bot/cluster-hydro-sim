#include "../src/moving_initialization.hpp"
#include "../src/moving_comparison.hpp"
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <cstring>
#include <chrono>

using S=MovingThreeFluidSim;
void check(bool condition,const char* message) {if(!condition)throw std::runtime_error(message);}
S setup(int n=32) {
  ThreeFluidSim hydro;hydro.initSolver(n);hydro.initPlummer(1,0.1,0.2,1,1);
  S s;s.param.epsilon=0.03;s.param.Deltat=1e-5;
  initializeMovingFromHydrostatic(s,hydro);return s;
}
Eigen::MatrixXd dense(const S& s) {
  int n=s.state.size();Eigen::MatrixXd a=Eigen::MatrixXd::Zero(n,n);
  for(int j=0;j<n;++j)for(int i=std::max(0,j-S::KU);i<=std::min(n-1,j+S::KL);++i)
    a(i,j)=s.bandMatrix()[S::KL+S::KU+i-j+S::LDAB*j];
  return a;
}
void bandSolve() {
  auto s=setup(12);s.param.c2={0.1,0.03,0.02};s.param.c1.fill(0.04);s.param.c4.fill(0.005);
  s.param.q=0.002;s.assembleStep();
  const auto a=dense(s);const auto rhs=Eigen::Map<const Eigen::VectorXd>(s.rightHandSide().data(),s.state.size()).eval();
  const Eigen::VectorXd expected=a.partialPivLu().solve(rhs);
  s.advanceAcceptedStep();const auto got=Eigen::Map<const Eigen::VectorXd>(s.lastIncrement().data(),s.state.size());
  check((a*got-rhs).norm()/(rhs.norm()+1e-300)<1e-10,"band residual");
  check((got-expected).norm()/(expected.norm()+1e-300)<1e-8,"band/dense disagreement");
}
void massLedger() {
  auto s=setup();std::array<double,3> mass;
  for(int f=0;f<3;++f)mass[f]=s.value(s.zones()-1,f,S::MASS);
  s.param.maxSteps=40;s.param.maxTime=1;s.param.c2={0.1,0.03,0.02};
  int calls=0;auto observer=[&](const S& x) {
    check(x.energy_ledger_error<1e-12,"energy/work ledger");
    ++calls;for(int f=0;f<3;++f)check(std::abs(x.value(x.zones()-1,f,S::MASS)+x.escaped_mass[f]-mass[f])<1e-12*mass[f],"mass ledger");
  };s.evolve(observer);check(calls==41,"observer calls");
  check(s.escaped_mass[0]>0&&s.escaped_heat[0]>0,"missing open fluxes");
}
void equilibriumAndTide() {
  auto s=setup(80);s.assembleStep();
  for(int i=0;i<s.zones()-1;++i)for(int f=0;f<3;++f) {
    check(std::abs(s.rightHandSide()[S::index(i,f,S::RHO)])<1e-20,"equilibrium mass residual");
    check(std::abs(s.rightHandSide()[S::index(i,f,S::U)])<1e-20,"equilibrium energy residual");
    check(std::abs(s.rightHandSide()[S::index(i,f,S::VEL)])<1e-14,"equilibrium momentum residual");
  }
  const auto original=s.rightHandSide();s.param.q=2;s.assembleStep();
  for(int i=0;i<s.zones();++i)for(int f=0;f<3;++f) {
    const double expected=s.Deltat*s.volume(i)*s.value(i,f,S::RHO)*s.param.q*s.radius(i)/s.param.epsilon;
    const double got=s.rightHandSide()[S::index(i,f,S::VEL)]-original[S::index(i,f,S::VEL)];
    check(std::abs(got-expected)<1e-12*expected,"tidal source or force-balance guard");
  }
}
void splitSymmetry() {
  auto s=setup(48);for(int i=0;i<s.zones();++i) {
    s.state[S::index(i,FD,S::RHO)]=s.value(i,FS,S::RHO);
    s.state[S::index(i,FD,S::U)]=s.value(i,FS,S::U);
  }
  auto faces=s.faces();auto state=s.state;s.initialize(faces,state,true);
  s.param.mass[FD]=s.param.mass[FS];s.param.c2[FD]=s.param.c2[FS]=0.17;
  s.param.c1[FS*3+FD]=s.param.c1[FD*3+FS]=0.325;
  s.param.maxSteps=25;auto obs=[](const S&){};s.evolve(obs);
  for(int i=0;i<s.zones();++i)for(int k:{S::RHO,S::U,S::MASS})
    check(std::abs(s.value(i,FS,k)/s.value(i,FD,k)-1)<1e-10,"equal-fluid symmetry");
}
void failEarly() {
  auto s=setup();s.param.epsilon=0;
  bool failed=false;try{s.validate();}catch(const std::exception&){failed=true;}
  check(failed,"epsilon validation");
  s=setup();s.Deltat=0;failed=false;
  try{s.advanceAcceptedStep();}catch(const std::exception&){failed=true;}
  check(failed&&s.step==0,"invalid timestep accepted");
  s=setup();s.state.push_back(1);failed=false;
  try{s.advanceAcceptedStep();}catch(const std::exception&){failed=true;}
  check(failed,"resized state reached band solve");
}
void dilutedReference() {
  auto s=setup();
  for(int i=s.zones()-4;i<s.zones();++i)for(int f=0;f<3;++f) {
    s.state[S::index(i,f,S::RHO)]*=1e-4;
    s.state[S::index(i,f,S::U)]*=0.01;
  }
  // This used to fail while constructing a face, before any solve: the
  // old equilibrium offset subtracted more density than remained locally.
  s.assembleStep();
  for(double x:s.bandMatrix())check(std::isfinite(x),"diluted reference coefficient");
}
void tidalEvolution() {
  auto s=setup();s.param.q=2;s.param.maxSteps=10;
  auto obs=[](const S&){};s.evolve(obs);
  check(s.value(0,FS,S::VEL)>0&&s.escaped_mass[0]>0,"outward q evolution");
}
void thermalAgreement() {
  ThreeFluidSim h;h.initSolver(32);h.initPlummer(1,0.1,0.2,1,1);
  h.param.md=0.1*h.param.ms;h.param.c1.fill(0.03);h.param.c4.fill(0.002);
  h.Deltat=1e-7;
  S s;s.param.epsilon=1;s.param.Deltat=h.Deltat;initializeMovingFromHydrostatic(s,h);
  h.solveConductionLAPACKE();s.advanceAcceptedStep();
  for(int i=3;i<20;++i)for(int f=0;f<3;++f)
    check(std::abs(s.value(i,f,S::U)/h.U[f][i]-1)<1e-11,"thermal stencil/legacy mismatch");
}
void parameterRoundTrip() {
  auto s=setup();s.param.q=.123;s.param.c1.fill(.314);s.param.mass={1,2,3};
  const std::string dir="output/check_moving_parameters/";
  std::filesystem::create_directories(dir);save_param_for_Mathematica(s.param,dir);
  MovingThreeFluidParam result{};std::ifstream file(dir+"param.dat",std::ios::binary);
  file.read(reinterpret_cast<char*>(&result),sizeof(result));
  check(file.gcount()==sizeof(result)&&std::memcmp(&result,&s.param,sizeof(result))==0,"parameter round trip");
}
void singleSplitAgreement() {
  std::array<S,2> sims;
  for(int k=0;k<2;++k) {
    ThreeFluidSim h;h.initSolver(40);h.initPlummer(k?0.5:1,1e-10,k?1:1e-10,1,1);
    h.param.md=h.param.ms;h.param.c2[FS]=h.param.c2[FD]=.17;
    sims[k].param.epsilon=1e-9;sims[k].param.Deltat=1e-5;sims[k].param.max_timestep=1e-5;
    sims[k].param.maxSteps=100;initializeMovingFromHydrostatic(sims[k],h);
    auto observe=[](const S&){};sims[k].evolve(observe);
  }
  for(int i=0;i<sims[0].zones();++i) {
    double total[2]{};
    for(int k=0;k<2;++k)for(int f=0;f<3;++f)total[k]+=sims[k].value(i,f,S::RHO);
    check(std::abs(total[0]/total[1]-1)<1e-8,"single/split total-density equivalence");
    check(std::abs(sims[0].value(i,FS,S::U)/sims[1].value(i,FS,S::U)-1)<1e-8,"single/split dispersion equivalence");
  }
}
void benchmark() {
  for(int n:{50,100,200,400,800}) {
    auto s=setup(n);s.param.maxSteps=50;s.param.max_timestep=s.param.Deltat;
    auto obs=[](const S&){};
    const auto start=std::chrono::steady_clock::now();s.evolve(obs);
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::cout<<"zones="<<n<<" seconds_per_step="<<seconds/s.step
      <<" band_bytes="<<s.bandMatrix().size()*sizeof(double)<<'\n';
  }
}
void exportFixture(const std::string& directory) {
  std::filesystem::create_directories(directory);
  constexpr int n=5;S s;s.param.epsilon=0.1;s.param.q=0.01;s.param.Deltat=0.01;
  s.param.c2.fill(0.07);s.param.c1.fill(0.003);s.param.c4.fill(0.002);
  s.param.mass={1,2,0.1};
  std::vector<double> faces{0,0.4,0.7,1.1,1.6,2.2},initial(12*n,0);
  for(int i=0;i<n;++i)for(int f=0;f<3;++f) {
    initial[S::index(i,f,S::RHO)]=1+0.1*f+0.04*i;
    initial[S::index(i,f,S::U)]=1+0.2*f+0.03*i;
    initial[S::index(i,f,S::VEL)]=0.2+0.03*f+0.01*i;
    if(i==n-1&&f>0) initial[S::index(i,f,S::VEL)]=(f==1?2:-4)*
      std::sqrt((10./9.)*initial[S::index(i,f,S::U)]/s.param.epsilon);
  }
  s.initialize(faces,initial,false);s.assembleStep();
  write_to_file(s.state,directory+"/state.dat");write_to_file(s.faces(),directory+"/faces.dat");
  const Eigen::MatrixXd a=dense(s);
  write_to_file(a,directory+"/matrix.dat");write_to_file(s.rightHandSide(),directory+"/rhs.dat");
}
void auditHydro() {
  std::cout<<std::setprecision(17);
  for(double dt:{1e-3,1e-4,1e-5,1e-6,1e-7,0.0}) {
    ThreeFluidSim s;MovingComparisonInit p;p.sample=2;initializeMovingComparison(s,p);
    s.param.c2.fill(0);s.Deltat=dt; // isolate local c1 exchange from conduction
    auto moments=[&]() {
      std::array<double,3> result{};
      for(int f=0;f<3;++f)for(int i=0;i<s.param.N;++i) {
        const double lo=i?s.R[f][i-1]:0,hi=s.R[f][i];
        const double mass=s.Rho[f][i]*(hi*hi*hi-lo*lo*lo)/3;
        result[0]+=mass;result[1]+=mass*s.U[f][i];
        result[2]+=mass*std::log(s.P[f][i]/std::pow(s.Rho[f][i],5./3.));
      }return result;
    };
    const auto initial=moments();
    const double u0=s.U[FS][0];
    s.solveConductionLAPACKE();const double heat_u=s.U[FS][0];
    const auto heated=moments();s.projectHydrostatic();const auto projected=moments();
    const double project_u=s.U[FS][0];s.realign();const auto aligned=moments();
    std::cout<<"dt="<<dt<<" align_dM="<<aligned[0]-projected[0]
      <<" align_dE="<<aligned[1]-projected[1]<<" align_dEntropy="<<aligned[2]-projected[2]
      <<" core_dU_heat_absolute="<<heat_u-u0;
    if(dt>0)std::cout<<" conduction_dE_dt="<<(heated[1]-initial[1])/dt
      <<" align_dM_dt="<<(aligned[0]-projected[0])/dt
      <<" align_dE_dt="<<(aligned[1]-projected[1])/dt
      <<" align_dEntropy_dt="<<(aligned[2]-projected[2])/dt
      <<" core_dU_heat="<<(heat_u-u0)/dt<<" core_dU_project="<<(project_u-heat_u)/dt
      <<" core_dU_align="<<(s.U[FS][0]-project_u)/dt;
    std::cout<<'\n';
  }
}
int main(int argc,char**argv) {try {
  if(argc==3&&std::string(argv[1])=="--export") {exportFixture(argv[2]);return 0;}
  if(argc==2&&std::string(argv[1])=="--audit-hydro") {auditHydro();return 0;}
  if(argc==2&&std::string(argv[1])=="--benchmark") {benchmark();return 0;}
  bandSolve();massLedger();equilibriumAndTide();splitSymmetry();failEarly();dilutedReference();tidalEvolution();
  thermalAgreement();parameterRoundTrip();singleSplitAgreement();
  std::cout<<"moving checks passed: band/dense, mass/heat export, equilibrium, q, split, fail-early\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
