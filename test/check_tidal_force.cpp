#include "../src/ab_reproduction.hpp"
#include <iostream>
#include <stdexcept>
#include <filesystem>
#include <limits>

void require(bool ok,const char* message) {if(!ok)throw std::runtime_error(message);}
double volume(const ThreeFluidSim& s,int f,int i) {
  return (std::pow(s.R[f][i],3)-(i?std::pow(s.R[f][i-1],3):0))/3;
}
// Independent residual, with the enclosed masses frozen as in assembly.
double residual(const ThreeFluidSim& s,int f,int i) {
  const double r=s.R[f][i], mass=s.Menc[FS][i]+s.Menc[FB][i]+s.Menc[FD][i];
  const double effective=mass-s.param.tidal_q*r*r*r;
  if(i==s.param.N-1)
    return -s.P[f][i]/(r-s.R[f][i-1])+s.Rho[f][i]*effective/(r*r);
  return 4*r*r*(s.P[f][i+1]-s.P[f][i])+effective*(s.Rho[f][i]+s.Rho[f][i+1])*
    (s.R[f][i+1]-(i?s.R[f][i-1]:0));
}
void perturb(ThreeFluidSim& s,const ThreeFluidSim& original,int f,int j,double delta) {
  s.R[f]=original.R[f];s.Rho[f]=original.Rho[f];s.P[f]=original.P[f];
  s.R[f][j]+=delta;
  for(int i=j;i<std::min(j+2,int(s.param.N));++i) {
    const double ratio=volume(original,f,i)/volume(s,f,i);
    s.Rho[f][i]*=ratio;s.P[f][i]*=std::pow(ratio,5./3);
  }
}
int main(){try {
  for(double q:{0.,1e-4}) {
    ThreeFluidSim s;s.param.tidal_q=q;s.initSolver(32);
    s.initPlummer(.5,1e-10,1,1,1,10);
    s.sanityCheck();s.assembleRelaxation(FS);
    const auto original=s;
    double max_error=0;
    for(int j=0;j<s.param.N;++j) {
      const double h=1e-5*(j?s.R[FS][j]-s.R[FS][j-1]:s.R[FS][0]);
      auto plus=original,minus=original;
      perturb(plus,original,FS,j,h);perturb(minus,original,FS,j,-h);
      for(int i=std::max(0,j-1);i<std::min(j+2,int(s.param.N));++i) {
        const double finite=(residual(plus,FS,i)-residual(minus,FS,i))/(2*h);
        const double analytic=i==j?original.hydroD[i]:(j<i?original.hydroDL[i-1]:original.hydroDU[i]);
        const double error=std::abs(finite-analytic)/std::max({std::abs(finite),std::abs(analytic),1e-30});
        max_error=std::max(max_error,error);
      }
    }
    require(max_error<1e-6,"tidal Jacobian finite-difference mismatch");
    s.solveRelaxationLAPACKE(FS);
    require(((s.R[FS]-original.R[FS]).array()/original.R[FS].array()).abs().maxCoeff()<1e-11,
            "tidal equilibrium changed under relaxation");
    for(int i=0;i<s.param.N;++i) {
      require(std::abs(s.Rho[FS][i]*volume(s,FS,i)/(original.Rho[FS][i]*volume(original,FS,i))-1)<1e-12,"shell mass");
      require(std::abs((s.P[FS][i]/std::pow(s.Rho[FS][i],5./3))/(original.P[FS][i]/std::pow(original.Rho[FS][i],5./3))-1)<1e-12,"entropy");
    }
    s.realign();s.sanityCheck();
    require(((s.U[FS]-original.U[FS]).array()/original.U[FS].array()).abs().maxCoeff()<1e-10,"tidal realignment fixed point");
    std::cout<<"q="<<q<<" max_jacobian_error="<<max_error<<" passed=1\n";
  }
  ThreeFluidSim uniform;uniform.initSolver(3);uniform.param.tidal_q=.002;
  for(int f=0;f<NF;++f){uniform.R[f]<<1,2,3;uniform.Rho[f].setConstant(1.0/NF);}
  require(std::abs(uniform.externalPotentialEnergy()+.002*std::pow(3.,5)/10)<1e-15,"external energy factor");
  bool rejected=false;
  try{ThreeFluidSim s;s.initSolver(32);s.param.tidal_q=.1;s.initPlummer(.5,1e-10,1,1,1,10);}
  catch(const std::runtime_error&){rejected=true;}
  require(rejected,"supercritical hydrostatic boundary accepted");
  for(double q:{-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
    ThreeFluidSim s;s.initSolver(32);s.initPlummer(.5,1e-10,1,1,1,10);s.param.tidal_q=q;
    rejected=false;try{s.validateEvolution();}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"invalid tidal coefficient accepted");
  }
  for(bool tide:{false,true}) {
    ThreeFluidSim s;s.param.N=80;ABInitParam init;auto obs=initializeAB(s,init,tide);
    require(obs.time_unit_myr>0&&obs.friction_time_myr>0,"AB unit conversion");
    const std::string dir=tide?"output/test_tidal_force/tidal/":"output/test_tidal_force/isolated/";
    std::filesystem::create_directories(dir);s.saveParams(dir);
    s.param.maxSteps=10;s.param.Deltat=1e-5;
    int count=0;auto observer=[&](const ThreeFluidSim&){++count;};s.evolve(observer);
    require(count==11,"AB evolve observation count");s.sanityCheck();
  }
  std::cout<<"tidal_force_check_passed=1\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
