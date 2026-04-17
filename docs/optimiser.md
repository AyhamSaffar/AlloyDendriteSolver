# Optimiser (WIP)

Each optimiser's objective is to find values for V and R that minimise the outputs of both f1 and f2 for any given
model.

This is done iteratively by picking values for V and R, calculating a tangent to the f surface at that point using
automatically generated gradients, and solving for the change in V and R at which those tangets would reach 0 f.

The Newton-Raphson scheme in two variables can be written in long form as:

$$ \begin{bmatrix} V_{n+1} \cr R_{n+1} \end{bmatrix} = \begin{bmatrix} V_n \cr R_n \end{bmatrix} - J^{-1}(V_n, R_n) \begin{bmatrix} f1(V_n, R_n) \cr f2(V_n, R_n) \end{bmatrix} $$

Where:
- n is the current step
- n+1 is the next step in the iteration
- $J$ is the Jacobian matrix:

$$ J= \begin{bmatrix} {\partial f_1}/{\partial V} & {\partial f_1}/{\partial R} \cr {\partial f_2}/{\partial V} & {\partial f_2}/{\partial R} \end{bmatrix} $$
