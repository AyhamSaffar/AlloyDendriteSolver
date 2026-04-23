# Models

Models are used to evaluate how accurate a given V and R pair are given C0, dT, and a given alloy. They are analytically
derived from theory given a set of physical assumptions.

When a given model is implemented, each equation given below is rearranged such that the right hand side equals zero.
This means that if all parameters are consistent, the right hand side of the first equation (f1) and the second equation
(f2) should evaluate to 0.

They must also be continuous so that they can be automatically differentiated at any given V and R pair.

The following alloy thermodynamic constants are used below:
- $L$ 	&nbsp; Latent heat of fusion - $J/kg$
- $c_p$ &nbsp; Specific heat capacity - $J/(kgK)$
- $m$ 	&nbsp; Equilibrium liquidus slope - $K/ wt.$%
- $k_0$	&nbsp; Partition coefficient - *unitless*
- $Γ$ 	&nbsp; Gibbs-Thomson coefficient - $Km$
- $D$ 	&nbsp; Solute diffusion coefficient - $m^2/s$
- $α$	&nbsp; Thermal diffusivity in the liquid - $m^2/s$
- $σ^*$	&nbsp; Stability constant - *unitless*
- $a_0$ &nbsp; Atomic spacing in pure metal - $m$
- $V_0$ &nbsp; Speed of sound in liquid - $m/s$
- $T_m$ &nbsp; Melting point of pure metal - $K$

- $R_0$ &nbsp; Gas constant - $J/(mol K)$

### The LGK model

[Lipton, J., Glicksman, M. E., & Kurz, W.](https://doi.org/10.1016/0025-5416(84)90199-X)

This equation holds up to moderate undercooling and when there is only a single nucleation event. The latter is common
in small liquid solder balls that don't have any available nucleants.

$$ ∆T = \frac{L}{c_p} Iv_t + mC_0 \{ 1 - \frac{1}{(1-(1-k_0 )Iv_c)} \} + \frac{2Γ}{R} $$

$$ R = \frac {Γ/σ^*} {\frac{P_t L}{c_p} -\frac{P_c m C_0 (1-k_0 )}{1-(1-k_0 ) Iv_c}} $$

Given the following

$P_t = \frac{VR}{2α}$ &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; - thermal Péclet number

$P_c = \frac{VR}{2D}$ &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; - solutal Péclet number

$Iv_t(P_t) = P_t e^{P_t} E_1(P_t)$ &emsp;&emsp;&emsp;&emsp;- solutal Ivantsov function

$Iv_c(P_c) = P_c e^{P_c} E_1(P_c)$ &emsp;&emsp;&emsp;&ensp;&nbsp;- thermal Ivantsov function

$E_1(x)$ &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;&emsp;&ensp;&nbsp;- first exponential integral of $x$

The first equation calculates the LGK dendrite undercooling. This equation takes into account the thermal,
solutal, and curvature undercooling. It uses dimensional analysis to solve for
solute and heat transport across an equilibrium solidification parabaloid dendrite and uses phase diagram
constants to calculate the drop in liquidus temperature ahead of the solidification front due to solute
enrichment

The second equation calculates the LGK stability criterion dendrite radius. A planar solidification front is modified by
adding a periodic pertubation, solving for the solute concentration field of the liquid side, and calculating the rate
of growth of a given wavelength with time. The smallest pertubation wavelength that does not grow with time can then be
used as the dendrite radius. This gives an expression in terms of temperature and solutal gradients, which can be
calculated for a parabaloid dendrite using the Ivantsov fields levered in the first equation.

### The LKT-BCT Model

[J. Lipton, W. Kurz, R. Trivedi](https://doi.org/10.1016/0001-6160(87)90174-X) - [W.J. Boettinger, S.R. Coriell and R. Trivedi (pages 13-31)](https://search.library.uq.edu.au/discovery/fulldisplay/alma991011497109703131/61UQ_INST:61UQ)

An extension of the LGK model that maintains accuracy at higher undercoolings and growth rates via fewer modelling
assumptions.

$$ ∆T = \frac{L}{c_p} Iv_t + mC_0 \{ 1 - \frac{m'/m}{(1-(1-k_0 )Iv_c)} \} + \frac{2Γ}{R} + \frac{V}{\mu} $$

$$ R = \frac {Γ/σ^*} {\frac{\xi_t P_t L}{c_p} -\frac{2 \xi_c P_c m' C_0 (1-k)}{1-(1-k) Iv_c}} $$

GIven the following:

$k = \frac{k_0 + (a_0V/D)}{1 + (a_0V/D)}$ &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; - velocity dependant partition coefficient

$m' = m(1 + \frac{k_0 - k[1-ln(k/k_0)]}{1-k_0})$ &emsp;&emsp;&emsp;&emsp; - velocity dependant liquidus slope

$\mu = \frac{LV_0}{R_0T_L^2}$ &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; - interfacial kinetic coefficient

$\xi_t = 1 - \frac{1}{\sqrt{1+1/(\sigma^*P_t^2)}}$ &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; - thermal stability function

$\xi_c = 1 + \frac{2k}{1-2k-\sqrt{1+1/(\sigma^*P_c^2)}}$ &emsp;&emsp;&emsp;&emsp;&emsp;&ensp; - solutal stability function

*explanation of the first equation*

*explanation of the second equation*
