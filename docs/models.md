# Models

TODO - list assumptions for each model. 

Models are used to evaluate how accurate a given V and R pair are given C0, dT, and a given alloy. They are analytically
derived from theory given a set of physical assumptions.

When a given model is implemented, each equation given below is rearranged such that the right hand side equals zero.
This means that if all parameters are consistent, the right hand side of the first equation (f1) and the second equation
(f2) should evaluate to 0.

They must also be continuous so that they can be automatically differentiated at any given V and R pair.

All alloy thermodynamic constants used below are documented in the [*alloys*](alloys.md) docs:

### The LGK model

[Lipton, J., Glicksman, M. E., & Kurz, W.](https://doi.org/10.1016/0025-5416(84)90199-X)

This equation holds up to moderate undercooling and when there is only a single nucleation event. The latter is common
in small liquid solder balls that don't have any available nucleants.

$$ ∆T = \frac{L}{c_p} Iv_t + mC_0 \left[ 1 - \frac{1}{1-(1-k_0 )Iv_c} \right] + \frac{2Γ}{R} $$

$$ R = \frac {Γ/σ^*} {\frac{P_t L}{c_p} -\frac{2 P_c m C_0 (1-k_0 )}{1-(1-k_0 ) Iv_c}} $$

Given the following

$P_t = \frac{VR}{2α}$ &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; - thermal Péclet number

$P_c = \frac{VR}{2D}$ &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; - solutal Péclet number

$Iv_t(P_t) = P_t e^{P_t} E_1(P_t)$ &emsp;&emsp;&emsp;&emsp;- solutal Ivantsov function

$Iv_c(P_c) = P_c e^{P_c} E_1(P_c)$ &emsp;&emsp;&emsp;&ensp;&nbsp;- thermal Ivantsov function

$E_1(x) = \int_u^\infty \frac{exp(-s)}{s} ds$ &emsp;&emsp;&ensp;&emsp;&ensp;&nbsp;- first exponential
integral of $x$

The first equation calculates the LGK dendrite undercooling. It quantifies how the liquid must be cooled below the
temperature of the solid to 1. drive thermal diffusion away from the solid that gives out heat as it solidifies, 2.
reach the lower melting temperature caused by a build up of solute just ahead of the solidification front, and 3.
overcome the energy barrier created by the surface energy of a high curvature dendrite tip. It uses dimensional analysis
to solve for solute and heat transport across an equilibrium solidification parabaloid dendrite. Phase diagram constants
are used to calculate the drop in liquidus temperature ahead of the solidification front due to solute enrichment.

The second equation calculates the LGK marginal stability criterion dendrite radius. A planar solidification front is
modified by adding a periodic pertubation. Too small and the curvature will drive the pertubation to shrink. Too large
and purtubation will grow by escaping the cold and solute rich solidification front. The dendrite radius is approximated
as the smallest pertubation that won't shrink. This gives an expression that is a function of the solute and temperature
field gradient, which can be calculated for a parabaloid dendrite using the same dimensional analysis as the first
equation.

Note the extra factor of 2 in the second term of the second equation's denominator. Lipton, Glicksman, & Kurz remove
this factor in their paper in order to coerce this equation into agreeing with a prior published result for the case
where there is zero thermal field gradient and the second equation only depends on the solutal field gradient. This
change is not otherwise justified and is ignored in future iterations of this model such as LKT-BCT.

*Assumptions* (unfinished)
- most constants do not vary with composition (E.G m, k0, r)

### The LKT-BCT Model

[J. Lipton, W. Kurz, R. Trivedi](https://doi.org/10.1016/0001-6160(87)90174-X) - [W.J. Boettinger, S.R. Coriell and R. 
Trivedi*](https://search.library.uq.edu.au/discovery/fulldisplay/alma991011497109703131/61UQ_INST:61UQ)

An extension of the LGK model that maintains accuracy at higher undercoolings and growth rates via fewer modelling
assumptions.

$$ ∆T = \frac{L}{c_p} Iv_t + mC_0 \left[ 1 - \frac{m'/m}{1-(1-k)Iv_c} \right] + \frac{2Γ}{R} + \frac{V}{\mu} $$

$$ R = \frac {Γ/σ^*} {\frac{\xi_t P_t L}{c_p} - \frac{2 P_c m C_0 (1-k) \xi_c}{1-(1-k) Iv_c}} $$

GIven the following:

$k = \frac{k_0 + (a_0V/D)}{1 + (a_0V/D) - (1-k_0)\frac{C_0}{100}}$ &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; - velocity
 dependant partition coefficient

$m' = m \left[1 + \frac{k_0 - k(1-ln(k/k_0))}{1-k_0} \right]$ &emsp;&emsp;&emsp;&emsp; - velocity dependant liquidus slope

$R_0$ = 8.314 $J/molK$ &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; - molar gas constant

$\mu = \frac{LV_0}{R_0T_m^2}$ &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; -
 interfacial kinetic coefficient

$\xi_t = 1 - \frac{1}{\sqrt{1+1/(\sigma^*P_t^2)}}$ &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; - thermal stability
 function

$\xi_c = 1 + \frac{2k}{1-2k-\sqrt{1+1/(\sigma^*P_c^2)}}$ &emsp;&emsp;&emsp;&emsp;&emsp;&ensp; - solutal stability
 function

Both terms modify their phase diagram constants ($m$ and $k$) to accomodate solute trapping at higher velocities. This
is when the solidification front moves too quickly to allow as many solute atoms to cross the solidification front into
the liquid. This means more solute gets frozen in the solid before crossing over to the liquid. At extreme velocities,
both solid and liquid have the same solute concentration ($C_0$), meaning $k \rightarrow$ 1 and both the liquidus and
solidus lines overlap with equal gradients.

The first equation has an added fourth undercooling term, which specifies how much further the liquid must be cooled to
overcome the kinetic energy barrier for adding extra liquid atoms onto the solid lattice. This is assumed to be
negligible at the lower velocities expected in the LGK model. 

The LGK second equation assumes small Peclet numbers, where V*R << 1, meaning the stability functions ≈ 1. The second
LKT-BCT equation however drops this low V assumption, meaning these terms must be included.

\* While Boettinger, Coriell, and Trivedi's origional paper is not openly published online, it is famously well written
and the basis for this implementation. For a full derivation of this model, it can be accessed by requesting pages 13-25
of the linked conference paper. Note that a $P_c$ is missing in the denominator of the paper's R equation, and the paper
assumes the dilute solute limit when deriving the velocity dependant $k$, meaning the $(1-k_0)C_0$ term is assumed to be
negligible. To better generalise to higher solute concentrations, this implementation keeps the term in.
