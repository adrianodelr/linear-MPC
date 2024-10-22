#include "linearMPC.h"

using namespace LMPC;

// SS model dimensions
const int n = 2;
const int m = 1;

// number of constraints 
const int pu = 1;       // number of control constraints 
const int px = 1;       // number of state constraints 

// horizon lenght 
const int N = 5;

// discrete SS model of a linear harmonic oscillator 
Matrix<n,n, float> A = {0.999,0.00999667,-0.199933,0.999};
Matrix<n,m, float> B = {0.0009998333444440899,0.19993333999968252};

// build model class object
auto model = LTIModel<n, m, float> (A, B); 

// limit vectors 
Matrix<m+m,1,float> ulimits = {10.0, 0.0/0.0};
Matrix<n+n,1,float> xlimits = {0.0/0.0, 15.0, 0.0/0.0, 0.0/0.0};

// build constraint object
auto constr = Constraints<n,m,pu,px,float>(ulimits,xlimits);

// Weight matrices 
Matrix<n,1,float> Q = {15.0,0.2};           // reference tracking weight  
Matrix<n,1,float> Qf = {15.0,0.2};          // reference tracking weight final
Matrix<m,1,float> R = {1.0};                // control weight

// desired setpoint/ state
Matrix<n,1,float> r_k = {10.0,0.0};

// build controller with constraints 
auto ctrl = lMPC<n,m,N,pu,px,float>(model, Q, Qf, R, r_k, constr); 

void setup() {
  Serial.begin(9600); 

  Matrix<n,1,float> x = {0,0}; // initial state  

  // simulate 50 time steps 
  int Nsim = 50; 
  for (int i = 0; i < Nsim; i++){
      auto u = ctrl.get_control(x);
      x = model.simulate_time_step(x, u);
      Serial.print("x: ");
      Serial.print(x);
      Serial.print("u: ");
      Serial.println(u);
  }
}

void loop() {
}
