# Linear Model Predictive Control

This is an easy to use Arduino library for controlling LTI systems in state-space form via Model Predictive Control (MPC). All matrix operations are handled by the [BasicLinearAlgebra library](https://github.com/tomstewart89/BasicLinearAlgebra) and the underlying optimization is done by the [Augmented Lagrangian Quadratic Program (QP) Solver](https://github.com/adrianodelr/ALQP-Solver). 

The MPC architecture is based on a condensed QP formulation, see e.g. [here, section 3a)](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=7ff6f36a6ff9a8461b11ea26bcc46a6db38443a6). This means, that the 
discrete-time dynamic equations (our state-space model), which are naturally linear equality constraints in the optimization, are substituted into the quadratic cost function. As a result,
future state variables are eliminated from the optimization search space, and the equality constraints are implicity satisfied. Instead of formulation the optimization problem in terms of absolute controls, this implementation uses **control increments**. The book of J.A. Rossiter [A First Course in Predictive Control](https://api.pageplace.de/preview/DT0400.9781351597166_A35143461/preview-9781351597166_A35143461.pdf) explains all these concepts very thoroughly for beginners. Detailed instructions for using the library are provided below.


## Usage
### Define a model 
The objective is to control a discrete-time state-space model of the form 

```math
    \mathbf{x}_{k+1} = \mathbf{A} \mathbf{x}_{k} + \mathbf{B} \mathbf{u}_{k}
```
with state vector $\mathbf{x}\in \mathbb{R}^{n}$, control vector $\mathbf{u}\in \mathbb{R}^{n}$. The user must provide the state matrix $\mathbf{A}\in \mathbb{R}^{n \times n}$ and the control matrix $\mathbf{B}\in \mathbb{R}^{n \times m}$:

```cpp
#include "linearMPC.h"

using namespace LMPC;

// SS model dimensions
const int n = 2;
const int m = 1;

// discrete SS model of a linear harmonic oscillator
Matrix<n,n,float> A = {0.999,0.00999667,-0.199933,0.999};
Matrix<n,m,float> B = {0.0009998333444440899,0.19993333999968252};

// build model class object
auto model = LTIModel<n, m, float> (A, B); 
```
### Define constraints (optional) 
Constraints can be used to reflect physical limits of the system. However, this increases computational time and can lead to a shortage of RAM. Limits can be set partially, on individual components of the state or the control vector.  

Upper and lower limits on controls are stacked together in a vector twice the size of the control vector, arranged in ascending order based on the corresponding vector elements, i.e. 

```math
    u_\text{limits} = \begin{bmatrix} u_{1max}, &  u_{2max}, & \dots & u_{m\; max}, &  u_{1min}, &  u_{2min}, &  \dots & u_{m \; min}\end{bmatrix}
```
The same applies to limits on the states:
```math
    x_\text{limits} = \begin{bmatrix} x_{1max}, &  x_{2max}, & \dots & x_{n\; max}, &  x_{1min}, &  x_{2min}, &  \dots & x_{n \; min}\end{bmatrix}
```
To set partial limits only on certain components of state or control vector, the corresponding value in the limits vector is set to `NaN`, which can be achieved by zero division, i.e. `0.0/0.0`. It is of utmost importance to correctly specify the number of constraints for each state and control vector.

In the example below, $u_{1max} = 10.0$  and $x_{2max} = 15.0$. There is no other limits, so $u_{1min} = x_{1max} = x_{1min} = x_{2min} = NaN$. Accordingly there is one constraint on the state and one on the controls.   

```cpp
// number of constraints 
const int pu = 1;       // number of control constraints 
const int px = 1;       // number of state constraints 

// limit vectors 
Matrix<m+m,1,float> ulimits = {10.0, 0.0/0.0};
Matrix<n+n,1,float> xlimits = {0.0/0.0, 15.0, 0.0/0.0, 0.0/0.0};

// build constraint object
auto constr = Constraints<n,m,pu,px,float>(ulimits,xlimits);
```
If there are no limits on either the states or the controls, you can pass a single vector to build the constraints object:
```cpp
// build constraint object with only control constraints
auto constr = Constraints<n,m,pu,0,float>(ulimits);

// OR

// build constraint object with only state constraints
auto constr = Constraints<n,m,0,px,float>(xlimits);
```
### Building a controller 
The unconstrained quadratic objective can be written using sum notation as
```math
    J = \mathbf{e}_{N}^T \mathbf{Q}_{f} \mathbf{e}_{N} + \sum_{k=1}^{N-1} \mathbf{e}_{k}^T \mathbf{Q} \mathbf{e}_{k} + \mathbf{u}_{k}^T \mathbf{R} \mathbf{u}_{k}

```

Here, $\mathbf{e}$<sub>$k$</sub> $=\mathbf{r}$<sub>$k$</sub> $-\mathbf{x}$<sub>$k$</sub> is the difference between the desired reference and the predicted state at time step $k$, and $\mathbf{u}_{k}$ the control vector applied to the system. The diagonal matrices $\mathbf{Q}$ and $\mathbf{R}$ determine how much emphasis is placed on reference tracking and control effort, respectively. The tracking error term at the last time step $N$ (the length of the prediction horizon) is taken out of the sum, as it is often given a higher value to emphasize the terminal state.  

For building a controller, only the diagonal terms of the weight matrices need to be specified. 
Additionally, the used must provide the desired reference $\mathbf{r}_k$. Currently, only **setpoint control** is supported, meaning $\mathbf{r}_k$ is constant.
```cpp
// horizon lenght 
const int hx = 5;

// Weight matrices 
Matrix<n,1,float> Q = {15.0,0.2};           // reference tracking weight  
Matrix<n,1,float> Qf = {15.0,0.2};          // reference tracking weight final
Matrix<m,1,float> R = {1.0};                // control weight

// desired setpoint/ state
Matrix<n,1,float> r_k = {10.0,0.0};

// build controller
// with constraints 
auto ctrl = lMPC<n,m,hx,pu,px,float>(model, Q, Qf, R, r_k, constr);  

// OR 
// without constraints
auto ctrl = lMPC<n,m,hx,0,0,float>(model, Q, Qf, R, r_k);  
```

The QP-solver parameters can be accessed via the `QPparams` struct which is closer described in the [ALQP-Solver documentation](https://github.com/adrianodelr/ALQP-Solver). The settings can be also useful for debugging purposes, as they allow to print available RAM during solving, along with some information about the solver's progress:

```cpp
// create params object
QPparams params;
params.debugging=true;
  
auto ctrl = lMPC<n,m,hx,pu,px,float>(model, Q, Qf, R, xref, constr, params);  
```

To use the controller in simulation, the `model` can be used: 

```cpp
Matrix<n,1,float> x = {0,0}; // initial state  

// simulate 15 time steps 
int Nsim = 15; 
for (int i = 0; i < Nsim; i++){
    auto u = ctrl.get_control(x);
    x = model.simulate_time_step(x, u);
}
```

