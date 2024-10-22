# Linear Model Predictive Control

This is an easy to use Arduino library for controlling LTI systems in state space form via Model Predictive Control (MPC). All matrix operations are handled by the [BasicLinerAlgebra library](https://github.com/tomstewart89/BasicLinearAlgebra) and the underlying optimization is done by the [Augmented Lagrangian Quadratic Program (QP) Solver](https://github.com/adrianodelr/ALQP-Solver). 

The MPC architecture is based on a [condensed QP formulation](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=7ff6f36a6ff9a8461b11ea26bcc46a6db38443a6). This essentially means, that the 
discrete time dynamic equations (here our state space model), which naturally are linear equality constraints
to be satisfied during the optimization, are substituted into the quadratic cost function. As a consequence,
future state variables are eliminated from the search space of the optimization, and the equality constraints
are implicity satisfied at the optimum. Instead formulation the optimization problem in terms of absolute controls, this implementation uses control increments. The book of J.A. Rossiter [A First Course in Predictive Control](https://api.pageplace.de/preview/DT0400.9781351597166_A35143461/preview-9781351597166_A35143461.pdf) explains all these concepts for beginners. However, for completeness, below the basics are layed out.   

## MPC formulation
$$
\begin{align}
\min_{\mathbf{x}} \quad & \frac{1}{2}\mathbf{x}^T\mathbf{Q}\mathbf{x} + \mathbf{q}^T\mathbf{x} \\ 
\mbox{s.t.}\quad &  \mathbf{A}\mathbf{x} -\mathbf{b} = \mathbf{0} \\ 
&  \mathbf{G}\mathbf{x} - \mathbf{h} \leq \mathbf{0} 
\end{align}
$$


### TEST FUTURE ERROR
$$ 
\underset{\rightarrow}{\mathbf{e}}_{k} 
$$
<!-- TEST FUTURE ERROR as command
$$
\newcommand{\efutr}[1] {\underset{\rightarrow}{\boldsymbol{e}}_{#1}}
$$

test of future vector format
$$
\efutr{k}  = 
\begin{bmatrix}
\boldsymbol{e}_{k+1} \\
\boldsymbol{e}_{k+2} \\
\vdots \\
\boldsymbol{e}_{k+n} \\
\end{bmatrix}
$$ -->


