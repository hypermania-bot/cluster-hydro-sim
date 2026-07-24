#pragma once

#include "Eigen/Dense"
#include <vector>
#include "param.hpp"
#include "io.hpp"



constexpr double PI = 3.14159265358979323846;
constexpr double C1 = 1.0 / sqrt(3.0 * PI);          // ≈0.326
constexpr double ALPHA = 2 * erf(sqrt(1.5)) - 2 * exp(-1.5) * sqrt(6.0 / PI);
constexpr double BETA_C = 0.45;
constexpr double C2 = sqrt(2.0) * pow(3.0, -2) * ALPHA * BETA_C; // ≈0.086
constexpr double GAMMA_LAMBDA = 0.4;
constexpr double GAMMA_COULOMB = 0.11;

// Fluid indices (avoid name clash with any variable)
constexpr int FS = 0;  // single stars
constexpr int FB = 1;  // binaries
constexpr int FD = 2;  // dark matter
constexpr int NF = 3;

struct ApproximateCentralDensityObserver {
  std::vector<double> rho_c_list;
  
  std::vector<double> t_list;
  std::vector<double> R_packed;
  std::vector<double> Rho_s_packed;
  std::vector<double> Rho_b_packed;
  std::vector<double> Rho_d_packed;
  std::vector<double> U_s_packed;
  std::vector<double> U_b_packed;
  std::vector<double> U_d_packed;

  ApproximateCentralDensityObserver(const std::vector<double> &rho_c_list_) : rho_c_list(rho_c_list_) {}
  
  void operator()(const Eigen::VectorXd &R, const std::vector<Eigen::VectorXd> &Rho, const std::vector<Eigen::VectorXd> &U, const double t) {
    const double rho_c = Rho[FS][0] + Rho[FB][0] + Rho[FD][0];
    const size_t current_idx = t_list.size();
    if(current_idx < rho_c_list.size() && rho_c >= rho_c_list[current_idx]) {
      t_list.push_back(t);
      R_packed.insert(R_packed.end(), R.begin(), R.end());
      Rho_s_packed.insert(Rho_s_packed.end(), Rho[FS].begin(), Rho[FS].end());
      Rho_b_packed.insert(Rho_b_packed.end(), Rho[FB].begin(), Rho[FB].end());
      Rho_d_packed.insert(Rho_d_packed.end(), Rho[FD].begin(), Rho[FD].end());
      U_s_packed.insert(U_s_packed.end(), U[FS].begin(), U[FS].end());
      U_b_packed.insert(U_b_packed.end(), U[FB].begin(), U[FB].end());
      U_d_packed.insert(U_d_packed.end(), U[FD].begin(), U[FD].end());
    }
  }
  
  void save(const std::string &dir) const {
    write_to_file(t_list, dir + "t.dat");
    write_to_file(R_packed, dir + "R.dat");
    write_to_file(Rho_s_packed, dir + "Rho_s.dat");
    write_to_file(Rho_b_packed, dir + "Rho_b.dat");
    write_to_file(Rho_d_packed, dir + "Rho_d.dat");
    write_to_file(U_s_packed, dir + "U_s.dat");
    write_to_file(U_b_packed, dir + "U_b.dat");
    write_to_file(U_d_packed, dir + "U_d.dat");
  }
};

struct ThreeFluidParam {
  // Evolution parameters
  long long int N;
  double ms;
  double mb;
  double md;
  std::array<double, NF> c2; // Conduction
  std::array<double, NF*NF> c1; // Dynamical heating
  std::array<double, NF*NF> c4;  // Binary heating

  // Numerical control parameters
  double Deltat;
  double StopDensity;
  double thres;
};

// ===================================================================
// Full 3-fluid simulation
// ===================================================================
class ThreeFluidSim {
public:
  // Evolution parameters
  long long int N = 200;
  double ms = 1.0;
  double mb = 2.0;
  double md = 1e-10;
  double c2[NF]; // Conduction
  double c1[NF][NF]; // Dynamical heating
  double c4[NF][NF]; // Binary heating

  // Numerical control parameters
  double Deltat = 1e-3;
  double StopDensity = 1e12;
  double thres = 1e-3;
  
  // Data
  double totalTime = 0.0;
  std::vector<Eigen::VectorXd> Rho, U, Menc, P, R;

  // Temporary variables
  std::vector<double> conductionAB;
  std::vector<double> conductionBSTD;
  std::vector<int> conductionIPIV;
  
  std::vector<Eigen::VectorXd> sqrtU;
  std::vector<Eigen::VectorXd> logRho;
  Eigen::VectorXd logR;

  Eigen::VectorXd hydroDL;
  Eigen::VectorXd hydroD;
  Eigen::VectorXd hydroDU;
  Eigen::VectorXd hydroB;

  Eigen::VectorXd newR;
  Eigen::VectorXd newRho;


  // Bulk output storage
  // std::vector<double> timeHistory;
  // std::vector<Eigen::VectorXd> RHistory;
  // std::vector<std::vector<Eigen::VectorXd>> RhoHistory, UHistory, MencHistory;

  std::string dir;

  // Initialize solver along with path for saving data
  ThreeFluidSim(const std::string &dir_);

  // Allocate and zero memory for the solver
  void initSolver(const int N);
  
  // Compute and assign coeffs c1[NF][NF], c2[NF] and c4[NF][NF]
  void initCoeffs();

  // Assign coeffs according to Yiming's paper
  void initCoeffsYiming();

  // Assign initial conditions using algorithm replicated from Yiming's paper
  void initPlummerYiming(const double xi1, const double xi2, const double zeta1, const double zeta2);
  // Assign initial conditions
  void initPlummer();

  void printParams() const;
  void printCoeffs() const;

  void saveParams() const;


  
  void updateEnclosedMass();

  // ----------------------------------------------------------------
  // 1. CONDUCTION + INTERACTION
  // ----------------------------------------------------------------
  void solveConductionLAPACKE();

  // ----------------------------------------------------------------
  // 2. RELAXATION (per fluid, hydrostatic + constant entropy)
  //    Exact copy of Mathematica SolveRelaxation + notes (26)-(28)
  // ----------------------------------------------------------------
  void solveRelaxationLAPACKE(const int f);
  
  // ----------------------------------------------------------------
  // 3. REALIGNMENT (common grid + interpolation)
  // ----------------------------------------------------------------
  void realign();

  // ----------------------------------------------------------------
  // 4. BINARY FORMATION (simple total-mass transfer, notes 1.4)
  // ----------------------------------------------------------------
  void applyBinaryFormation(const double dt);

  // ----------------------------------------------------------------
  // 5. MAIN EVOLUTION (with adaptive time steps)
  // ----------------------------------------------------------------
  // void evolve(const int maxSteps);

  // template<typename Observer>
  void evolve(const int maxSteps, ApproximateCentralDensityObserver &observer);

  
};

