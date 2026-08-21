#include "three_fluid.hpp"
#include "io.hpp"
#include "utility.hpp"


#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <limits>
#include <lapacke.h>


using namespace Eigen;
using namespace std;


// ===================================================================
// Full 3-fluid simulation
// ===================================================================

ThreeFluidSim::ThreeFluidSim() {}


void ThreeFluidSim::initSolver(const int N) {
  this->N = N;
  for(int f = 0; f < NF; ++f){
    R[f].resize(N);
    Rho[f].resize(N);
    U[f].resize(N);
    Menc[f].resize(N);
    P[f].resize(N);
    sqrtU[f].resize(N);
    logRho[f].resize(N);
  }
  
  logR = VectorXd::Zero(N);

  // Init conduction solver
  {
    const int n_eq = NF * N;       // total number of unknowns
    const int kl = NF;                 // number of subdiagonals (3)
    const int ku = NF;                 // number of superdiagonals (3)
    const int ldab = 2 * kl + ku + 1;      // leading dimension of banded storage (7)
    conductionAB = vector<double>(ldab * n_eq, 0.0);
    conductionBSTD = vector<double>(n_eq);
    conductionIPIV = vector<int>(n_eq);
  }

  // Init relaxation solver
  {
    hydroDL = VectorXd::Zero(N-1);
    hydroD = VectorXd::Zero(N);
    hydroDU = VectorXd::Zero(N-1);
    hydroB = VectorXd::Zero(N);
  }

  // Init realign
  newR = VectorXd::Zero(N);
  newRho = VectorXd::Zero(N);
}


void ThreeFluidSim::initCoeffs() {
  // c' coefficients (exact from your notes eqs 51-57)
  constexpr double beta[NF] = {1.0, 1.0, 1.0};
  constexpr double MtotSys = 1e9;
  double lnLsd = log(GAMMA_LAMBDA * 2 * MtotSys / (ms + md));
  double mi[NF] = {ms, mb, md};
  for (int i = 0; i < NF; ++i) {
    double lnLi = log(GAMMA_LAMBDA * 2 * MtotSys / (mi[i] + mi[i]));
    c2[i] = 2.0 * C2 * beta[i] * mi[i] / ms * lnLi / lnLsd;
  }
  for (int i = 0; i < NF; ++i) {
    for (int j = 0; j < NF; ++j) {
      if(i != j){
	double lnLij = log(GAMMA_LAMBDA * 2 * MtotSys / (mi[i] + mi[j]));
	c1[i][j] = lnLij / lnLsd * C1;
      } else {
	c1[i][j] = 0;
      }
    }
  }

  c4[FS][FS] = 0;
  c4[FS][FB] = 0.123679 / lnLsd * ms / (mb + ms);
  c4[FS][FD] = 0;
    
  c4[FB][FS] = 0.123679 / lnLsd * ms / mb * ms / (mb + ms);
  c4[FB][FB] = 5.0 / 16.0 * sqrt(6.0 * PI) / lnLsd * mb / ms;
  c4[FB][FD] = 0.123679 / lnLsd * (md * md / mb / ms) * (md / (mb + md));

  c4[FD][FS] = 0;
  c4[FD][FB] = 0.123679 / lnLsd * md / ms * md / (mb + md);
  c4[FD][FD] = 0;

  binary_formation = 0;
  tidal_cutoff = 0;
}


void ThreeFluidSim::initCoeffsYiming() {
  // THIS BLOCK IS MANUALLY ADDED TO SET COEFFICIENTS THE SAME AS THAT IN YIMING'S CODE
  // Coefficients set up match Yiming's two fluid code
      
  // No binary heating
  c4[FS][FS] = 0;
  c4[FS][FB] = 0;
  c4[FS][FD] = 0;
    
  c4[FB][FS] = 0;
  c4[FB][FB] = 0;
  c4[FB][FD] = 0;

  c4[FD][FS] = 0;
  c4[FD][FB] = 0;
  c4[FD][FD] = 0;

  c2[FS] = 0.1720743853447534;
  c2[FB] = 0;
  c2[FD] = 0;

  c1[FS][FD] = c1[FD][FS] = 0.32573500793527993;
  c1[FS][FB] = c1[FB][FS] = 0;
  c1[FB][FD] = c1[FD][FB] = 0;
  c1[FS][FS] = c1[FB][FB] = c1[FD][FD] = 0;

  binary_formation = 0;
  tidal_cutoff = 0;
}


void ThreeFluidSim::printCoeffs() const {
  cout << "Initiated coefficients to:" << endl;
  cout << endl;
  for(int i = 0; i < NF; ++i){
    cout << "c2[" << i << "] = " << c2[i] << endl;
  }
  cout << endl;
  
  for(int i = 0; i < NF; ++i){
    for(int j = 0; j < NF; ++j){
	cout << "c4[" << i << "][" << j << "] = " << c4[i][j] << endl;
    }
  }
  cout << endl;
  for(int i = 0; i < NF; ++i){
    for(int j = 0; j < NF; ++j){
	cout << "c1[" << i << "][" << j << "] = " << c1[i][j] << endl;
    }
  }
}

void ThreeFluidSim::printParams() const {
  cout << "Initiated coefficients to:" << endl;
  cout << "N = " << N << endl;
  cout << "ms = " << ms << endl;
  cout << "mb = " << mb << endl;
  cout << "md = " << md << endl;

  cout << "binary_formation = " << binary_formation << endl;
  cout << "tidal_cutoff = " << tidal_cutoff << endl;
}


void ThreeFluidSim::initPlummerYiming(const double rho0, const double xi1, const double xi2, const double zeta1, const double zeta2){
  constexpr double r0 = 1.0;
  //constexpr double rho0 = 1.0;
  
  // Replicating Yiming's code
  const int layer = 150;
  const int extralayer = 10;
  N = layer + extralayer;
  for(int f = 0; f < NF; ++f) {
    R[f].resize(N);
    Rho[f].resize(N);
    U[f].resize(N);
    P[f].resize(N);
    Menc[f].resize(N);
    R[f].array() = pow(10.0, VectorXd::LinSpaced(N, -2.0, 3.0).array());
  }
  double rss = r0, rsb = zeta1 * r0, rsd = zeta2 * r0;
  double rhoss = rho0, rhosb = xi1 * rho0, rhosd = xi2 * rho0;
  for (int i = 0; i < N; ++i) {
    double r2 = (i==0) ? pow(R[0][0] / 2, 2) : pow((R[0][i-1] + R[0][i]) / 2, 2);
    // Assuming zeta_1, zeta_2 != 1
    Rho[FS][i] = rhoss / pow(1.0 + r2 / (rss * rss), 2.5);
    Rho[FB][i] = rhosb / pow(1.0 + r2 / (rsb * rsb), 2.5);
    Rho[FD][i] = rhosd / pow(1.0 + r2 / (rsd * rsd), 2.5);
  }
  // Rho[FS][N-1] = Rho[FB][N-1] = Rho[FD][N-1] = 0.0;
  
  updateEnclosedMass();



  // P[FS][N-1] = 0;
  // P[FB][N-1] = 0;
  // P[FD][N-1] = 0;

  // To be replaced with more accurate analytic solutions
  P[FS][N-1] = (1.0 + xi1 + xi2) / pow(1.0 + pow((R[0][N-2] + R[0][N-1]) / 2, 2), 3) / 18.0;
  P[FB][N-1] = xi1 * P[FS][N-1];
  P[FD][N-1] = xi2 * P[FS][N-1];
  for(int i = N-2; i >= 0; --i){
    double r_avg = (i==0) ? (R[0][0] / 2) : ((R[0][i-1] + R[0][i]) / 2);
    double r_next = (R[0][i] + R[0][i+1]) / 2;
    double M_tot = Menc[FS][i] + Menc[FB][i] + Menc[FD][i];
    P[FS][i] = P[FS][i+1] + M_tot / pow(R[0][i], 2) * ((Rho[FS][i] + Rho[FS][i+1]) / 2.0) * (r_next - r_avg);
    P[FB][i] = P[FB][i+1] + M_tot / pow(R[0][i], 2) * ((Rho[FB][i] + Rho[FB][i+1]) / 2.0) * (r_next - r_avg);
    P[FD][i] = P[FD][i+1] + M_tot / pow(R[0][i], 2) * ((Rho[FD][i] + Rho[FD][i+1]) / 2.0) * (r_next - r_avg);
  }

  for(int f = 0; f < NF; ++f){
    U[f].array() = 1.5 * P[f].array() / Rho[f].array();
  }

  // Replicating Yiming's code
  N = layer;
  VectorXd temp = VectorXd::Zero(N);
  for(int f = 0; f < NF; ++f) {
    temp = R[f].head(N);
    R[f] = temp;
    temp = Rho[f].head(N);
    Rho[f] = temp;
    temp = U[f].head(N);
    U[f] = temp;
    temp = P[f].head(N);
    P[f] = temp;
    temp = Menc[f].head(N);
    Menc[f] = temp;
  }

  updateEnclosedMass();

  // Reconstruct pressure after truncation so that zone N-1 is an active
  // finite-mass shell satisfying the same one-sided equation as relaxation.
  for(int f = 0; f < NF; ++f) {
    const double outer_width = R[f][N-1] - R[f][N-2];
    const double outer_mass = Menc[FS][N-1] + Menc[FB][N-1] + Menc[FD][N-1];
    P[f][N-1] = outer_mass * Rho[f][N-1] * outer_width / pow(R[f][N-1], 2);
  }
  for(int i = N-2; i >= 0; --i) {
    for(int f = 0; f < NF; ++f) {
      const double radial_span = (i == 0) ? R[f][1] : R[f][i+1] - R[f][i-1];
      const double enclosed_mass = Menc[FS][i] + Menc[FB][i] + Menc[FD][i];
      P[f][i] = P[f][i+1] + enclosed_mass * (Rho[f][i] + Rho[f][i+1])
        * radial_span / (4.0 * pow(R[f][i], 2));
    }
  }
  for(int f = 0; f < NF; ++f) {
    U[f].array() = 1.5 * P[f].array() / Rho[f].array();
  }
}
  
void ThreeFluidSim::initPlummer(const double rho0, const double xi1, const double xi2, const double zeta1, const double zeta2) {
  constexpr double r0 = 1.0;
  
  for(int f = 0; f < NF; ++f) {
    R[f].array() = pow(10.0, VectorXd::LinSpaced(N, -2.0, 3.0).array());
  }
    
  double rss = r0, rsb = zeta1 * r0, rsd = zeta2 * r0;
  double rhoss = rho0, rhosb = xi1 * rho0, rhosd = xi2 * rho0;
  for (int i = 0; i < N; ++i) {
    double r2 = (i==0) ? pow(R[0][0] / 2, 2) : pow((R[0][i-1] + R[0][i]) / 2, 2);
    Rho[FS][i] = rhoss / pow(1.0 + r2 / (rss * rss), 2.5);
    Rho[FB][i] = rhosb / pow(1.0 + r2 / (rsb * rsb), 2.5);
    Rho[FD][i] = rhosd / pow(1.0 + r2 / (rsd * rsd), 2.5);
    U[FS][i] = rhoss * pow(rss, 3) / pow(1.0 + r2 / (rss * rss), 0.5) / 12;
    U[FB][i] = rhosb * pow(rsb, 3) / pow(1.0 + r2 / (rsb * rsb), 0.5) / 12;
    U[FD][i] = rhosd * pow(rsd, 3) / pow(1.0 + r2 / (rsd * rsd), 0.5) / 12;
    // P[FS][i] = (2.0 / 3.0) * U[FS][i] * Rho[FS][i];
    // P[FB][i] = (2.0 / 3.0) * U[FB][i] * Rho[FB][i];
    // P[FD][i] = (2.0 / 3.0) * U[FD][i] * Rho[FD][i];
  }
  // cout << "Initiated profile to:" << endl;
  // cout << "Rho[FS] = " << Rho[FS].transpose() << endl;
  // cout << "U[FS] = " << U[FS].transpose() << endl;
  // cout << "P[FS] = " << P[FS].transpose() << endl;
  // cout << endl;

  updateEnclosedMass();
  for(int f = 0; f < NF; ++f) {
    const double outer_width = R[f][N-1] - R[f][N-2];
    const double outer_mass = Menc[FS][N-1] + Menc[FB][N-1] + Menc[FD][N-1];
    P[f][N-1] = outer_mass * Rho[f][N-1] * outer_width / pow(R[f][N-1], 2);
  }
  for(int i = N-2; i >= 0; --i) {
    for(int f = 0; f < NF; ++f) {
      const double radial_span = (i == 0) ? R[f][1] : R[f][i+1] - R[f][i-1];
      const double enclosed_mass = Menc[FS][i] + Menc[FB][i] + Menc[FD][i];
      P[f][i] = P[f][i+1] + enclosed_mass * (Rho[f][i] + Rho[f][i+1])
        * radial_span / (4.0 * pow(R[f][i], 2));
    }
  }
  for(int f = 0; f < NF; ++f) {
    U[f].array() = 1.5 * P[f].array() / Rho[f].array();
  }
}

void ThreeFluidSim::updateEnclosedMass() {
  for (int f = 0; f < NF; ++f) {
    Menc[f][0] = Rho[f][0] * pow(R[f][0], 3) / 3.0;
    for (int i = 1; i < N; ++i) {
      Menc[f][i] = Menc[f][i-1] + Rho[f][i] * (pow(R[f][i], 3) - pow(R[f][i-1], 3)) / 3.0;
    }
  }
}

// ----------------------------------------------------------------
// CONDUCTION + INTERACTION
// THIS IS THE VERSION OF CODE solveConduction WE ACTUALLY USE
// ----------------------------------------------------------------
void ThreeFluidSim::solveConductionLAPACKE() {
  // const int n_rad = N;               // number of radial points per fluid
  // const int nf = NF;                 // number of fluids (3)
  const int n_eq = NF * N;       // total number of unknowns
  const int kl = NF;                 // number of subdiagonals (3)
  const int ku = NF;                 // number of superdiagonals (3)
  const int ldab = 2 * kl + ku + 1;      // leading dimension of banded storage (7)

  auto &conductionB = conductionBSTD;
    
  for (int f = 0; f < NF; ++f) {
    sqrtU[f] = U[f].array().sqrt().matrix();
    logRho[f] = Rho[f].array().log().matrix();
  }
  logR = R[FS].array().log().matrix();

  // ---------- Build the banded matrix and RHS ----------
  fill(conductionAB.begin(), conductionAB.end(), 0.0);

  // Helper lambda to set an entry in banded storage (column‑major)
  auto set_band = [&](const int row, const int col, const double val)->void{
    int k = ku + kl + row - col;          // row index inside the band
    // conductionAB[k + col * ldab] = val;
    // // -ku <= row - col <= kl
    // // kl <= ku + kl + row - col <= 2 * kl + ku < ldab
    if (k >= kl && k < ldab) {
      conductionAB[k + col * ldab] = val;
    } else {
      cout << "(solveConductionLAPACKE) Incorrect banded matrix index!" << endl;
      exit(0);
    }
  };


    
  // ----- Innermost boundary (i = 0) -----
  // Enforces dU/dr = 0 boundary condition
  for (int f = 0; f < NF; ++f) {
    int row = f;                     // i=0, fluid f
    set_band(row, row, 1.0);         // diagonal
    set_band(row, row + NF, -1.0);   // coupling to i=1, same fluid
    conductionB[row] = 0.0;
  }
  // ----- Outermost boundary (i = N-1) -----
  // Currently simply fixes U[N-1], equivalent to an external heat source
  // Heat loss is typically negligible due to low temperature at outer boundary
  // Should be changed to a L=0 boundary condition
  for (int f = 0; f < NF; ++f) {
    int row = (N - 1) * NF + f;
    set_band(row, row, 1.0);
    conductionB[row] = U[f][N - 1];
  }

  double mi[3] = {ms, mb, md};   // masses (FS, FB, FD)

  // ----- Bulk zones: i = 1 ... N-2 -----
  for (int i = 1; i < N - 1; ++i) {
    double dlogR = logR[i] - logR[i-1];
    double R2dlogR2 = dlogR * dlogR * R[FS][i-1] * R[FS][i];

    for (int f = 0; f < NF; ++f) {
      int row = i * NF + f;

      // Coefficients for same fluid (tridiagonal part)
      double coefL = c2[f] * Deltat * (2*logR[i] - 2*logR[i-1] - 4 - logRho[f][i-1] + logRho[f][i+1]) / (8 * sqrtU[f][i-1]) / R2dlogR2;
      double coefR = c2[f] * Deltat * (2*logR[i-1] - 2*logR[i] - 4 + logRho[f][i-1] - logRho[f][i+1]) / (8 * sqrtU[f][i+1]) / R2dlogR2;
      double coefC = 1.0 + c2[f] * Deltat / sqrtU[f][i] / R2dlogR2;

      // Cross‑fluid terms (same i)
      for (int f2 = 0; f2 < NF; ++f2) {
	coefC += Deltat * c1[f][f2] * (mi[f] / ms) * Rho[f2][i] / pow(U[f][i] + U[f2][i], 1.5)
	  + Deltat * c4[f][f2] * Rho[f2][i] / (2 * pow(U[f][i], 1.5));
      }

      // Set band entries
      set_band(row, row - NF, coefL);   // previous i, same f
      set_band(row, row, coefC);        // diagonal
      set_band(row, row + NF, coefR);   // next i, same f

      // Cross terms (different fluids, same i)
      for (int f2 = 0; f2 < NF; ++f2) {
	if (f2 == f) continue;
	int col = i * NF + f2;
	double cross = - Deltat * c1[f][f2] * (mi[f2] / ms) * Rho[f2][i] / pow(U[f][i] + U[f2][i], 1.5)
	  + Deltat * c4[f][f2] * Rho[f2][i] / (2 * pow(U[f2][i], 1.5));
	set_band(row, col, cross);
      }

      {
	// set RHS
	double val = U[f][i]
	  + (c2[f] * Deltat) / (8.0 * R2dlogR2)
	  * ( -8 * sqrtU[f][i]
	      + sqrtU[f][i-1] * (2*logR[i-1] - 2*logR[i] + 4 + logRho[f][i-1] - logRho[f][i+1])
	      + sqrtU[f][i+1] * (2*logR[i] - 2*logR[i-1] + 4 - logRho[f][i-1] + logRho[f][i+1]) );
	for (int f2 = 0; f2 < NF; ++f2) {
	  val += 1.5 * c4[f][f2] * Deltat * Rho[f2][i] / (2 * pow(U[f][i], 1.5));
	}
	conductionB[i * NF + f] = val;
      }
    }
  }

  // ----- Build right‑hand side (explicit diffusion) -----
  // Reorder RHS into the same i‑major, f‑minor layout.
  // for (int f = 0; f < NF; ++f) {

  //   // bulk i = 1 ... N-3
  //   for (int i = 1; i < N - 2; ++i) {
  // 	double dlogR = logR[i] - logR[i-1];
  // 	double R2dlogR2 = dlogR * dlogR * R[FS][i-1] * R[FS][i];
  // 	double val = U[f][i]
  // 	  + (c2[f] * Deltat) / (8.0 * R2dlogR2)
  // 	  * ( -8 * sqrtU[f][i]
  // 	      + sqrtU[f][i-1] * (2*logR[i-1] - 2*logR[i] + 4 + logRho[f][i-1] - logRho[f][i+1])
  // 	      + sqrtU[f][i+1] * (2*logR[i] - 2*logR[i-1] + 4 - logRho[f][i-1] + logRho[f][i+1]) );
  // 	for (int f2 = 0; f2 < NF; ++f2) {
  // 	  val += 1.5 * c4[f][f2] * Deltat * Rho[f2][i] / (2 * pow(U[f][i], 1.5));
  // 	}
  // 	conductionB[i * NF + f] = val;
  //   }

  //   // i = N-2
  //   {
  // 	int i = N - 2;
  // 	double dlogR = logR[i] - logR[i-1];
  // 	double R2dlogR2 = dlogR * dlogR * R[FS][i-1] * R[FS][i];
  // 	double val = U[f][i]
  // 	  + (c2[f] * Deltat) / (8.0 * R2dlogR2)
  // 	  * ( -8 * sqrtU[f][i]
  // 	      + sqrtU[f][i-1] * (2*logR[i-1] - 2*logR[i] + 4 + 2*logRho[f][i-1] - 2*logRho[f][i])
  // 	      + sqrtU[f][i+1] * (2*logR[i] - 2*logR[i-1] + 4 - 2*logRho[f][i-1] + 2*logRho[f][i]) );
  // 	for (int f2 = 0; f2 < NF; ++f2) {
  // 	  val += 1.5 * c4[f][f2] * Deltat * Rho[f2][i] / (2 * pow(U[f][i], 1.5));
  // 	}
  // 	conductionB[i * NF + f] = val;
  //   }

  // }

  // ---------- Solve using LAPACKE_dgbsv ----------
  int info = LAPACKE_dgbsv(LAPACK_COL_MAJOR,   // storage order
			   n_eq,               // number of equations
			   kl, ku,             // bandwidth
			   1,                  // number of RHS
			   conductionAB.data(), ldab,    // banded matrix
			   conductionIPIV.data(),        // pivot array (output)
			   conductionB.data(), n_eq);    // RHS (overwritten with solution)

  if (info != 0) {
    std::cerr << "Conduction failed!" << info << std::endl;
    std::cerr << "LAPACKE_dgbsv failed with error code " << info << std::endl;
    exit(0);
  }

  // ---------- Unpack solution back to U[f][i] ----------
  // Solution is now in b[] with ordering: i major, f minor.
  for (int i = 0; i < N; ++i) {
    for (int f = 0; f < NF; ++f) {
      double val = conductionB[i * NF + f];
      U[f][i] = val;
    }
  }

  // Update P
  for(int f = 0; f < NF; ++f) {
    P[f] = (2.0 / 3.0) * (Rho[f].array() * U[f].array()).matrix();
  }
  
}


// ----------------------------------------------------------------
// 2. RELAXATION (per fluid, hydrostatic + constant entropy)
//    Exact copy of Mathematica SolveRelaxation + notes (26)-(28)
//    WE USE THIS CODE
// ----------------------------------------------------------------
void ThreeFluidSim::solveRelaxationLAPACKE(const int f) {
  const int nrhs = 1;
  const int ldb = N;

  // Innermost (i=0)
  double Mtot_i = Menc[FS][0] + Menc[FB][0] + Menc[FD][0];
  hydroD[0] = 8 * R[f][0] * (P[f][1] - P[f][0]) + 20 * pow(R[f][0], 4) *
    (P[f][0] / pow(R[f][0], 3) + P[f][1] / (pow(R[f][1], 3) - pow(R[f][0], 3))) +
    3 * Mtot_i * (R[f][1] - 0) * pow(R[f][0], 2) *
    (-Rho[f][0] / pow(R[f][0], 3) + Rho[f][1] / (pow(R[f][1], 3) - pow(R[f][0], 3)));
  hydroDU[0] = Mtot_i * (Rho[f][0] + Rho[f][1]) - 20 * pow(R[f][0], 2) * P[f][1] * pow(R[f][1], 2) / (pow(R[f][1], 3) - pow(R[f][0], 3)) -
    3 * Mtot_i * Rho[f][1] * (R[f][1]) * pow(R[f][1], 2) / (pow(R[f][1], 3) - pow(R[f][0], 3));
  hydroB[0] = -4 * pow(R[f][0], 2) * (P[f][1] - P[f][0]) - Mtot_i * (Rho[f][0] + Rho[f][1]) * R[f][1];

  // Bulk zones
  for (int i = 1; i < N - 1; ++i) {
    Mtot_i = Menc[FS][i] + Menc[FB][i] + Menc[FD][i];
    hydroDL[i-1] = -Mtot_i * (Rho[f][i] + Rho[f][i+1]) -
      20 * pow(R[f][i], 2) * P[f][i] * pow(R[f][i-1], 2) / (pow(R[f][i], 3) - pow(R[f][i-1], 3)) +
      3 * Mtot_i * (R[f][i+1] - R[f][i-1]) * Rho[f][i] * pow(R[f][i-1], 2) / (pow(R[f][i], 3) - pow(R[f][i-1], 3));
    hydroD[i] = 8 * R[f][i] * (P[f][i+1] - P[f][i]) +
      20 * pow(R[f][i], 4) * (P[f][i] / (pow(R[f][i], 3) - pow(R[f][i-1], 3)) +
			      P[f][i+1] / (pow(R[f][i+1], 3) - pow(R[f][i], 3))) +
      3 * Mtot_i * (R[f][i+1] - R[f][i-1]) * pow(R[f][i], 2) *
      (-Rho[f][i] / (pow(R[f][i], 3) - pow(R[f][i-1], 3)) + Rho[f][i+1] / (pow(R[f][i+1], 3) - pow(R[f][i], 3)));
    hydroDU[i] = Mtot_i * (Rho[f][i] + Rho[f][i+1]) -
      20 * pow(R[f][i], 2) * P[f][i+1] * pow(R[f][i+1], 2) / (pow(R[f][i+1], 3) - pow(R[f][i], 3)) -
      3 * Mtot_i * Rho[f][i+1] * (R[f][i+1] - R[f][i-1]) * pow(R[f][i+1], 2) / (pow(R[f][i+1], 3) - pow(R[f][i], 3));
    hydroB[i] = -4 * pow(R[f][i], 2) * (P[f][i+1] - P[f][i]) - Mtot_i * (Rho[f][i] + Rho[f][i+1]) * (R[f][i+1] - R[f][i-1]);
  }

  // Outermost one-sided equation:
  // -P[N-1] / (R[N-1] - R[N-2])
  //   + Mtot[N-1] * Rho[N-1] / R[N-1]^2 = 0.
  // The Jacobian preserves the outer shell's mass and specific entropy.
  const int last = N - 1;
  const double inner_radius = R[f][last-1];
  const double outer_radius = R[f][last];
  const double outer_density = Rho[f][last];
  const double outer_pressure = P[f][last];
  const double outer_mass = Menc[FS][last] + Menc[FB][last] + Menc[FD][last];
  const double outer_volume = pow(outer_radius, 3) - pow(inner_radius, 3);
  const double outer_width = outer_radius - inner_radius;

  hydroDL[last-1] =
    -5 * outer_pressure * pow(inner_radius, 2) / (outer_volume * outer_width)
    - outer_pressure / pow(outer_width, 2)
    + 3 * outer_mass * outer_density * pow(inner_radius, 2)
      / (outer_volume * pow(outer_radius, 2));
  hydroD[last] =
    5 * outer_pressure * pow(outer_radius, 2) / (outer_volume * outer_width)
    + outer_pressure / pow(outer_width, 2)
    - 3 * outer_mass * outer_density / outer_volume
    - 2 * outer_mass * outer_density / pow(outer_radius, 3);
  hydroB[last] = outer_pressure / outer_width
    - outer_mass * outer_density / pow(outer_radius, 2);

  int info = LAPACKE_dgtsv(LAPACK_COL_MAJOR, N, nrhs,
			   hydroDL.data(),
			   hydroD.data(),
			   hydroDU.data(),
			   hydroB.data(),
			   ldb);

  if(info != 0){
    cout << "Relaxation failed!" << endl;
    cout << "hydroDL = " << hydroDL.transpose() << endl;
    cout << "hydroD = " << hydroD.transpose() << endl;
    cout << "hydroDU = " << hydroDU.transpose() << endl;
    cout << "hydroB = " << hydroB.transpose() << endl;
    exit(0);
  }

  bool R_ordered = true;
  double R_i = R[f][0] + hydroB[0];
  for(int i = 1; i < N; ++i){
    double next_R_i = R[f][i] + hydroB[i];
    if(R_i >= next_R_i){
      R_ordered = false;
      break;
    }
    R_i = next_R_i;
  }

  if(!R_ordered){
    cout << "Relaxation: unordered R[f] after relaxation!" << endl;
    exit(0);
  }

  // Version that exactly preserves s_i and m_i
  {
    double s_i = P[f][0] / pow(Rho[f][0], 5.0/3.0);
    Rho[f][0] = Rho[f][0] * pow(R[f][0], 3) / pow(R[f][0] + hydroB[0], 3);
    P[f][0] = s_i * pow(Rho[f][0], 5.0/3.0);
    
    U[f][0] = 1.5 * P[f][0] / Rho[f][0];
    for (int i = 1; i < N; ++i) {
      s_i = P[f][i] / pow(Rho[f][i], 5.0/3.0);
      Rho[f][i] = Rho[f][i] * (pow(R[f][i], 3) - pow(R[f][i-1], 3)) / (pow(R[f][i] + hydroB[i], 3) - pow(R[f][i-1] + hydroB[i-1], 3));
      P[f][i] = s_i * pow(Rho[f][i], 5.0/3.0);
      U[f][i] = 1.5 * P[f][i] / Rho[f][i];
    }
  }

  // Version that use linear approximation
  // {
  //   Rho[f][0] += -3 * Rho[f][0] * hydroB[0] / R[f][0];
  //   P[f][0] += -5 * P[f][0] * hydroB[0] / R[f][0];
    
  //   U[f][0] = 1.5 * P[f][0] / Rho[f][0];
  //   for (int i = 1; i < n; ++i) {
  // 	double ratio = (R[f][i] * R[f][i] * hydroB[i] - R[f][i-1] * R[f][i-1] * hydroB[i-1]) / (pow(R[f][i], 3) - pow(R[f][i-1], 3));
  // 	Rho[f][i] += -3 * Rho[f][i] * ratio;
  // 	P[f][i] += -5 * P[f][i] * ratio;
	
  // 	U[f][i] = 1.5 * P[f][i] / Rho[f][i];
  //   }
  //   Rho[f][N-1] = 0;
  //   P[f][N-1] = 0;
  //   U[f][N-1] = U[f][N-2];
  // }

  R[f] += hydroB;
}
  
// ----------------------------------------------------------------
// 3. REALIGNMENT (common grid + interpolation)
// ----------------------------------------------------------------
void ThreeFluidSim::realign() {
  double rMin = min({R[FS][0], R[FB][0], R[FD][0]});
  double rMax = max({R[FS][N-1], R[FB][N-1], R[FD][N-1]});
    
  newR = rMin * pow(rMax / rMin, VectorXd::LinSpaced(N, 0.0, 1.0).array()).matrix();
		      
    
  for (int f = 0; f < NF; ++f) {
    // Align density (log-log interp)
    for (int i = 0; i < N; ++i) {
      auto it = lower_bound(R[f].begin(), R[f].end(), newR[i]);
      int idx = it - R[f].begin();
      if(idx == 0) {
	// Assuming dRho/dr = 0 in first Lagrangian zone
	newRho[i] = Rho[f][0]; 
      } else if(idx < N) {
	// Logarithmic interpolation between R[f][idx-1] and R[f][idx]
	double t = (log(newR[i]) - log(R[f][idx-1])) / (log(R[f][idx]) - log(R[f][idx-1]));
	newRho[i] = exp(log(Rho[f][idx-1]) * (1-t) + log(Rho[f][idx]) * t);
	// cout << "f,i,idx,newR,t,newRho,Rho = " << f << "," << i << "," << idx << "," << newR[i] << "," << t << "," << newRho[i] << "," << Rho[f][i] << endl;
      } else {
	// Extrapolate from the two outermost active shells.
	double t = (log(newR[i]) - log(R[f][N-2])) / (log(R[f][N-1]) - log(R[f][N-2]));
	newRho[i] = exp(log(Rho[f][N-2]) * (1-t) + log(Rho[f][N-1]) * t);
	// cout << "f,i,idx,newR,t,newRho,Rho = " << f << "," << i << "," << idx << "," << newR[i] << "," << t << "," << newRho[i] << "," << Rho[f][i] << endl;
      }
    }
    Rho[f] = newRho;
      
    // Recompute Menc
    Menc[f][0] = newRho[0] * pow(newR[0], 3) / 3.0;
    for (int i = 1; i < N; ++i) {
      Menc[f][i] = Menc[f][i-1] + newRho[i] * (pow(newR[i], 3) - pow(newR[i-1], 3)) / 3.0;
    }
      
  }
  for(int f = 0; f < NF; ++f){
    R[f] = newR;
  }



  // Recompute U from hydrostatic (AlignU style)
  for (int f = 0; f < NF; ++f) {
    const double outer_width = newR[N-1] - newR[N-2];
    const double outer_mass = Menc[FS][N-1] + Menc[FB][N-1] + Menc[FD][N-1];
    P[f][N-1] = outer_mass * Rho[f][N-1] * outer_width / pow(newR[N-1], 2);
    for (int i = N-2; i >= 1; --i) {
      //double avgRho = (Rho[f][i] + Rho[f][i+1]) / 2.0;
      P[f][i] = P[f][i+1] + (Menc[FS][i] + Menc[FB][i] + Menc[FD][i]) * (Rho[f][i] + Rho[f][i+1]) * (newR[i+1] - newR[i-1]) / (4.0 * pow(newR[i], 2));
      // newU[i] = 1.5 * pShell / Rho[f][i];
    }
    P[f][0] = P[f][1] + (Menc[FS][0] + Menc[FB][0] + Menc[FD][0]) * (Rho[f][0] + Rho[f][1]) * newR[1] / (4.0 * pow(newR[0], 2));

    U[f] = 1.5 * (P[f].array() / Rho[f].array()).matrix();
  }
}

// ----------------------------------------------------------------
// 4. BINARY FORMATION
// ----------------------------------------------------------------
void ThreeFluidSim::applyBinaryFormation() {
  
  if(binary_formation == 1) {
    // (simple total-mass transfer, notes 1.4)
    using namespace std::numbers;
    
    // Calculate total binary formation rate
    double dM_dt = 0;
    for(int i = 0; i < N; ++i){
      const double MtotSys = 1e9;
      const double lnLsd = log(GAMMA_LAMBDA * 2 * MtotSys / (ms + md));

      // dimensionless formation rate from tidal capture
      // dn_dt is the time derivative of binary number density
      // For now assum ms = M_sun, Rs = R_sun, v_m = 10km/s
      // double dn_dt_tc = 131513.0 * pow(Rho[FS][i], 2) / lnLsd;
      // double dn_dt_tc = pow(Rho[FS][i], 2) / lnLsd;
      double dn_dt_tc = 0;

      // dimensionless formation rate from 3-body interaction
      double dn_dt_3b = 1e-8 * 1.65256 * ms * pow(Rho[FS][i], 3) / (lnLsd * pow(U[FS][i], 4.5));
      // double dn_dt_3b = 0;
      double vol = (i == 0) ? (pow(R[FS][i], 3) / 3.0) : ((pow(R[FS][i], 3) - pow(R[FS][i-1], 3)) / 3.0);
      dM_dt += ms * (dn_dt_tc + dn_dt_3b) * vol;
    }

    // Scale Rho so that mass is transferred from single star to binary
    double M_FS_ratio = (Menc[FS][N-1] - dM_dt * Deltat) / Menc[FS][N-1];
    double M_FB_ratio = (Menc[FB][N-1] + dM_dt * Deltat) / Menc[FB][N-1];
    for(int i = 0; i < N; ++i){
      Rho[FS][i] *= M_FS_ratio;
      Rho[FB][i] *= M_FB_ratio;
    }
    updateEnclosedMass();

    for(int f = 0; f < NF; ++f) {
      P[f].array() = (2.0 / 3.0) * (Rho[f].array() * U[f].array());
    }
    
    std::cout << "M_FS_ratio, M_FB_ratio, Menc[FS], Menc[FB] = "
	      << M_FS_ratio << "," << M_FB_ratio << "," << Menc[FS][N-1] << "," << Menc[FB][N-1] << std::endl;
  } else if(binary_formation == 2) {
    // Transfer by each Lagrangian zone
    // To preserve energy conservation at each binary formation step, we require
    // \begin{align}
    // & \rho_s u_s + \rho_b u_b = \rho_s' u_s' + \rho_b' u_b' \nonumber \\
    // & u_s = u_s'
    // \end{align}
    // The first equation fixes total heat energy at each Lagrangian zone.
    // The second equation states that binary formation doesn't change the internal energy of the single star specie.
    
    for(int i = 0; i < N; ++i){
      const double MtotSys = 1e9;
      const double lnLsd = log(GAMMA_LAMBDA * 2 * MtotSys / (ms + md));

      // dimensionless formation rate from tidal capture
      // dn_dt is the time derivative of binary number density
      // For now assum ms = M_sun, Rs = R_sun, v_m = 10km/s
      // double dn_dt_tc = 131513.0 * pow(Rho[FS][i], 2) / lnLsd;
      // double dn_dt_tc = pow(Rho[FS][i], 2) / lnLsd;
      double dn_dt_tc = 0;

      // dimensionless formation rate from 3-body interaction
      double dn_dt_3b = 1e-9 * 1.65256 * ms * pow(Rho[FS][i], 3) / (lnLsd * pow(U[FS][i], 4.5));
      // double dn_dt_3b = 0;
      // double vol = (i == 0) ? (pow(R[FS][i], 3) / 3.0) : ((pow(R[FS][i], 3) - pow(R[FS][i-1], 3)) / 3.0);
      double dRho = mb * (dn_dt_tc + dn_dt_3b) * Deltat;

      double RhoB_new = (Rho[FS][i] * U[FS][i] + Rho[FB][i] * U[FB][i] - (Rho[FS][i] - dRho) * U[FS][i]) / (Rho[FB][i] + dRho);
      
      Rho[FS][i] -= dRho;
      Rho[FB][i] += dRho;
      U[FB][i] = RhoB_new;
    }
    updateEnclosedMass();

    for(int f = 0; f < NF; ++f) {
      P[f].array() = (2.0 / 3.0) * (Rho[f].array() * U[f].array());
    }
  }
  
}


void ThreeFluidSim::applyTidalCutoff() {
  //const double tidal_radius = min({R[FS][N-1], R[FB][N-1], R[FD][N-1]});
  // const double tidal_radius = 1e1;
  
  for(int f = 0; f < NF; ++f) {
    // if(tidal_radius < R[f][0]) {
    //   Rho[f][0] *= pow(tidal_radius / R[f][0], 3);
    // }
    for(int i = 1; i < N; ++i) {
      // This doesn't work because of numerical instability
      // if(tidal_radius < R[f][i-1]) {
      // 	Rho[f][i] = 1e-10;
      // } else if(tidal_radius < R[f][i]) {
      // 	double old_vol = pow(R[f][i], 3) - pow(R[f][i-1], 3);
      // 	double new_vol = pow(tidal_radius, 3) - pow(R[f][i-1], 3);
      // 	Rho[f][i] *= new_vol / old_vol;
      // }

      // This is equivalent to a - const * rho term.
      // d rho / d t = - factor * rho
      const double factor = tidal_cutoff_factor;
      if(tidal_radius < R[f][i-1]){
	Rho[f][i] *= 1.0 - factor * Deltat;
      } else if(tidal_radius < R[f][i]) {
	const double t = (tidal_radius - R[f][i-1]) / (R[f][i] - R[f][i-1]);
	Rho[f][i] *= 1.0 - (1.0 - t) * factor * Deltat;
      }
    }
    P[f].array() = (2.0 / 3.0) * U[f].array() * Rho[f].array();
  }
  updateEnclosedMass();
}


void ThreeFluidSim::saveParams(const std::string &dir) const {
  ThreeFluidParam param;
  param.N = N;
  param.ms = ms;
  param.mb = mb;
  param.md = md;
  std::memcpy(param.c2.data(), c2, sizeof(c2));
  std::memcpy(param.c1.data(), c1, sizeof(c1));
  std::memcpy(param.c4.data(), c4, sizeof(c4));
  param.binary_formation = binary_formation;
  param.tidal_cutoff = tidal_cutoff;
  param.tidal_cutoff_factor = tidal_cutoff_factor;
  param.tidal_radius = tidal_radius;
  param.Deltat = Deltat;
  param.StopDensity = StopDensity;
  param.thres = thres;

  save_param_for_Mathematica(param, dir);
}
