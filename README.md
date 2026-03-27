# Alloy Solidification Solver

We are interested in developing code that runs two models of dendrite growth into an undercooled melt, generating data linking four variables over a wide range of growth undercoolings and velocities:

    ΔT	Growth undercooling
    V	Tip velocity
    R	Tip radius
    C0	Bulk alloy composition

The LGK model considers dendrite growth into an undercooled melt at velocities where the solid-liquid interface remains at local equilibrium.

The LGK-BCT model considers dendrite growth into an undercooled melt at velocities where rapid solidification effects are significant, including a velocity dependent partition coefficient and velocity dependent liquidus slope.  

### The LGK model (1984)
[Lipton, J., Glicksman, M. E., & Kurz, W. (1984). Dendritic growth into undercooled alloy metals. Materials Science and Engineering, 65(1), 57–63.](https://doi.org/10.1016/0025-5416(84)90199-X)

$$ ∆T = \frac{L}{c_p} Iv_t + mC_0 \{ 1 - \frac{1}{(1-(1-k_0 )Iv_c)} \} + \frac{2Γ}{R} $$

$$ R = \frac {\frac{Γ}{σ^*}} {\frac{P_t L}{c_p} -\frac{P_c m C_0 (1-k_0 )}{1-(1-k_0 ) Iv_c}} $$

With the following identities

$P_t = \frac{VR}{2α}$ 

$P_c = \frac{VR}{2D}$

$Iv_t(P_t) = P_t e^{P_t} E_1(P_t)$

$Iv_c(P_c) = P_c e^{P_c} E_1(P_c)$

Where
- $P_t$       &nbsp; is the thermal Péclet number
- $P_c$       &nbsp; is the solutal Péclet number
- $Iv_c(P_c)$ &nbsp; is the solutal Ivantsov function
- $Iv_t(P_t)$ &nbsp; is the thermal Ivantsov function
- $E_1(x)$    &nbsp; is the first exponential integral of variable $x$

Within these equations, there are eight constant parameters, defined as follows:
- $L$ 	&nbsp; Latent heat of fusion - $J/kg$
- $c_p$ &nbsp; Specific heat capacity - $J/(kgK)$
- $m$ 	&nbsp; Equilibrium liquidus slope - $K/ wt.$%
- $k_0$	&nbsp; Partition coefficient - *unitless*
- $Γ$ 	&nbsp; Gibbs-Thomson coefficient - $Km$
- $D$ 	&nbsp; Solute diffusion coefficient - $m^2/s$
- $α$	&nbsp; Thermal diffusivity in the liquid - $m^2/s$
- $σ^*$	&nbsp; Stability constant - *unitless*

For Sn-Ag alloys, these parameters are
- $L$ = 	  61,810.62 $J/kg$
- $c_p$ =     249 $J/(kgK)$
- $m$ = 	  −3.14 $K/ wt.$%
- $k_0$ =	  0.0191
- $Γ$ = 	  8.54 * $10^8$ $Km$
- $D$ = 	  1.82 * $10^{–9}$ $m^2/s$
- $α$ =	      1.5 * $10^{–5}$ $m^2/s$
- $σ^*$ =	  $1/(4π^2)$

There are four variables, defined as follows:

- $C_0$ Bulk alloy composition - $wt.$%
- $∆T$ Undercooling - $K$
- $V$ Velocity - $m/s$
- $R$ Dendrite tip radius - $m$

Fixing two of these variables, the other two can be solved iteratively using a method such as a two-dimensional Newton scheme, also known as the Newton-Raphson method in two variables.

In the codes we have written in the past, the user defines fixed values for $C_0$ and $∆T$ and the values of the 8 constant parameters.  The code then uses the Newton-Raphson method in two variables to converge on the unique values of V and R that satisfy these two functions. For this, we need to get the LGK equations into a form equal to zero and define them as functions f1 and f2:

$$ f_1(V, R) = \frac{L}{c_p} Iv_t + mC_0 \{ 1 - \frac{1}{(1-(1-k_0 )Iv_c)} \} + \frac{2Γ}{R} -∆T = 0 $$

$$ f_2(V, R) = \frac {\frac{Γ}{σ^*}} {\frac{P_t L}{c_p} -\frac{P_c m C_0 (1-k_0 )}{1-(1-k_0 ) Iv_c}} - R = 0 $$

The Newton-Raphson scheme in two variables can be written in long form as:

$$ \begin{bmatrix} V_{n+1} \cr R_{n+1} \end{bmatrix} = \begin{bmatrix} V_n \cr R_n \end{bmatrix} - J^{-1}(V_n, R_n) \begin{bmatrix} f1(V_n, R_n) \cr f2(V_n, R_n) \end{bmatrix} $$

Where:
- n is the current step
- n+1 is the next step in the iteration
- $J$ is the Jacobian matrix:

$$ J= \begin{bmatrix} {\partial f_1}/{\partial V} & {\partial f_1}/{\partial R} \cr {\partial f_2}/{\partial V} & {\partial f_2}/{\partial R} \end{bmatrix} $$

To perform the Newton-Raphson method, we need a good initial guess or approximate solution to the values of $V_0$ and $R_0$, i.e. $V_n$ and $R_n$ when n=0.

This can be taken from the approximate analytical solution in Eq. 8.91 and 8.92 of the book ‘Solidification’ by Dantzig & Rappaz (1st Ed):

$$ R = 6.64π^2 Γ(-m(1-k_0 ))^{0.25} \frac{C_0^{0.25}}{ΔT^{1.25}} $$

$$ V = \frac{D} {5.51π^2 (-m(1-k_0 ))^{1.5}Γ} \frac{ΔT^{2.5}}{C_0^{1.5}} $$

This is only a good approximation when:
1. The undercooling is small
1. The dendrites are solutal

### The LKT-BCT model (1987 – 1988)

- Lipton J, Kurz W, Trivedi R. Rapid Dendrite Growth in Undercooled Alloys. Acta Metall. 1987;35(4):957–64.
- W.J. Boettinger, S.R. Coriell, R. Trivedi, Application of dendritic growth theory to the interpretation of rapid solidification microstructures, Rapid Solidif. Process. Princ. Technol. IV (1988) 13.

The equations are written out in simple form in the following papers:
- Appendix to: Sun, S., Li, A., Cheng, C., & Gourlay, C. M. (2025). Effects of Ag and melt undercooling on the microstructure of Sn–Ag solder balls. Journal of Materials Science: Materials in Electronics, 36(942). https://doi.org/10.1007/s10854-025-14979-6
- Rodriguez, J. E., Kreischer, C., Volkmann, T., & Matson, D. M. (2017). Solidification velocity of undercooled Fe-Co alloys. Acta Materialia, 122, 431–437. https://doi.org/10.1016/j.actamat.2016.09.047


### Existing Graphical Datasets To Compare To

For the LGK model:
- Fig. 3, 4, 5 in Lipton J, Glicksman ME, Kurz W. Dendritic growth into undercooled alloy metals. Mater Sci Eng. 1984 Jul;65(1):57–63.
- Fig. 14 in Boettinger WJ, Bendersky LA, Early JG. An analysis of the microstructure of rapidly solidified Al-8 wt pct Fe powder. Metall Trans A. 1986 May;17: 781–90.

For the LKT-BCT model:
- Fig. 2 in Boettinger WJ, Coriell SR, Trivedi R. Application of Dendritic Growth Theory to the Interpretation of Rapid Solidification Microstructures. In: Rapid Solidification Processing: Principles and Technologies IV. 1988. p. 13–25.
- Herlach, D. M., Eckler, K., Karma, A., & Schwarz, M. (2001). Grain refinement through fragmentation of dendrites in undercooled melts. Materials Science and Engineering: A, 304–306(1–2), 20–25. https://doi.org/10.1016/S0921-5093(00)01553-7
    - Specifically, you want to use your implementation of the LKT-BCT model to calculate R versus ΔT using the constants in their Table 1.  And then put your values for R versus ΔT into their Eq. 5 to plot out Δt_bu versus ΔT and see if you get their Figure 3.  Note in their Eq. 5 that, when they write R(ΔT)^3, I think they mean (R as a function of temperature) cubed…

### To Do
- Currently getting unreliable convergence
    - Maybe play with adaptive / clamped updates
    - Hess / more sophisticated optimizer
    - Maybe try looking up different techniques
    - Plot V vs R for fixed C0 & dT to see if multiple solutions (might explain weird converged results)
    - Find convergence intervals for initials guesses in V vs dT with lines for approx, LGK, and eventually LKT-BCT. This will be helpful in checking if fancier optimisers / update schemes are worthwhile in terms on interval of convergence.
- Attempt to reproduce some LGK plots
- Add T dependant LGK model where diffusivity, liquidus temperature, partition coefficient, and liquidus slopes are fit polynomials. *must figure out how best to implement this for stability criterion part of LGK*.
- Implement LKT-BCT model and reproduce published plots
- Fix README by moving model explanations into docs folder and adopting the summary, method, installation, usage, support, roadmap, and acknowledgment structure.
