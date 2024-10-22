# Linear Model Predictive Control

This is an easy to use Arduino library for controlling LTI systems in state space form via Model Predictive Control (MPC). All matrix operations are handled by the [BasicLinerAlgebra library](https://github.com/tomstewart89/BasicLinearAlgebra) and the underlying optimization is done by the [Augmented Lagrangian Quadratic Program (QP) Solver](https://github.com/adrianodelr/ALQP-Solver). 

The MPC architecture is based on a [condensed QP formulation](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=7ff6f36a6ff9a8461b11ea26bcc46a6db38443a6). This essentially means, that the 
discrete time dynamic equations (here our state space model), which naturally are linear equality constraints
to be satisfied during the optimization, are substituted into the quadratic cost function. As a consequence,
future state variables are eliminated from the search space of the optimization, and the equality constraints
are implicity satisfied at the optimum. Instead formulation the optimization problem in terms of absolute controls, this implementation uses control increments. The book of J.A. Rossiter [A First Course in Predictive Control](https://api.pageplace.de/preview/DT0400.9781351597166_A35143461/preview-9781351597166_A35143461.pdf) explains all these concepts for beginners. However, for completeness, the most important mathematics are layed out at the end of the readme.   



## Usage
### Define model 
The goal is to control a discrete time state space model of the form 

```math
    \mathbf{x}_{k+1} = \mathbf{A} \mathbf{x}_{k} + \mathbf{B} \mathbf{u}_{k}
```
with state vector $\mathbf{x}\in \mathbb{R}^{n}$, control vector $\mathbf{u}\in \mathbb{R}^{n}$. The user has to provide the state matrix $\mathbf{A}\in \mathbb{R}^{n \times n}$ and the control matrix $\mathbf{B}\in \mathbb{R}^{n \times m}$:

```cpp
#include "linearMPC.h"

using namespace LMPC;

// SS model dimensions
const int n = 2;
const int m = 1;

// discrete SS model of linear harmonic oscillator (for more info see end of readme)
Matrix<n,n,float> A = {0.999,0.00999667,-0.199933,0.999};
Matrix<n,m,float> B = {0.0009998333444440899,0.19993333999968252};

// build model class object
auto model = LTIModel<n, m, float> (A, B); 
```
### Define constraints (optional) 
Constraints can be used to reflect physical limits of the system. However, this increases computational time and can lead to a shortage of RAM. Limits can be set partically, on individual components of the state or the control vector. Upper and lower limits are stacked together 

```cpp
// number of constraints 
const int pu = 1;
const int px = 1;

// limits on 
Matrix<m+m,1,float> ulimits = {10.0,0.0/0.0};
Matrix<n+n,1,float> xlimits = {0.0/0.0,15,0.0/0.0,0.0/0.0};

// build constraint class object
  auto constr = Constraints<n,m,pu,px,float>(ulimits,xlimits);
```


## MPC formulation
$$
\begin{align}
\min_{\mathbf{x}} \quad & \frac{1}{2}\mathbf{x}^T\mathbf{Q}\mathbf{x} + \mathbf{q}^T\mathbf{x} \\ 
\mbox{s.t.}\quad &  \mathbf{A}\mathbf{x} -\mathbf{b} = \mathbf{0} \\ 
&  \mathbf{G}\mathbf{x} - \mathbf{h} \leq \mathbf{0} 
\end{align}
$$


### TEST FUTURE ERROR
Some rendering tests  

$$ 
\underset{\rightarrow k}{\mathbf{e}} 
$$  

The prediction matrix  

$$
\underset{\rightarrow k}{\mathbf{e}} =
$$

new line   
```math
  \begin{bmatrix} 
    X \\ 
    Y 
  \end{bmatrix}
```

newer line  

```math
\underbrace{\begin{bmatrix}
\boldsymbol{r}_{k+1} \\
\boldsymbol{r}_{k+2} \\
\vdots \\
\boldsymbol{r}_{k+n_y}
\end{bmatrix}}_\text{$\boldsymbol{\bar{A}}$}
```

newest line  

$$\begin{bmatrix}  \boldsymbol{r}_{k+1} &  \boldsymbol{r}_{k+1} \\\  \boldsymbol{r}_{k+1} &  \boldsymbol{r}_{k+1} \end{bmatrix}$$  


$$\begin{bmatrix} \boldsymbol{r}_{k+1} \\\ \boldsymbol{r}_{k+2} \\\ \vdots \\\ \boldsymbol{r}_{k+n_y} \end{bmatrix}$$


<!-- % -         
% \underbrace{\begin{bmatrix}
% \mathbf{A} \\
% \mathbf{A}^2 \\
% \vdots \\
% \mathbf{A}^{n_y}
% \end{bmatrix}}_\text{$\boldsymbol{\bar{A}}$}
% \mathbf{x}_{k} -
% \underbrace{\begin{bmatrix}
% \mathbf{B} & \mathbf{0} & \mathbf{0} & \cdots   \\
% \mathbf{A}\mathbf{B} & \mathbf{B} & \mathbf{0} & \cdots\\
% \vdots & \vdots & \vdots & \ddots\\
% \mathbf{A}^{n_y-1}\mathbf{B} & \mathbf{A}^{n_y-2}\mathbf{B} & \mathbf{A}^{n_y-3}\mathbf{B} & \vdots\\
% \end{bmatrix}}_\text{$\boldsymbol{\bar{B}}$}
% \underbrace{\begin{bmatrix}
% \boldsymbol{r}_{k} \\
% \boldsymbol{r}_{k+1} \\
% \vdots \\
% \boldsymbol{r}_{k+n_y-1}
% \end{bmatrix}}_\text{$\underset{\rightarrow k}{\boldsymbol{r}}$} -->