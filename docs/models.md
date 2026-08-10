# Models

Models are used to evaluate how accurate a given V and R pair are given C0, dT, and a given alloy. They are analytically
derived from theory given a set of physical assumptions. All models assume a single nucleation event, which is common
in small liquid solder balls that don't have any available nucleants.

When a given model is implemented, each equation given below is rearranged such that the right hand side equals zero.
This means that if all parameters are consistent, the right hand side of the first equation (f1) and the second equation
(f2) should evaluate to 0.

They must also be continuous so that they can be automatically differentiated at any given V and R pair.

All alloy thermodynamic constants used below are documented in the [*alloys*](alloys.md) docs:

### The LGK model

[Lipton, J., Glicksman, M. E., & Kurz, W.](https://doi.org/10.1016/0025-5416(84)90199-X)

This equation holds up to moderate undercooling for alloys with linear phase diagrams. 

$$ ∆T = \frac{L}{c_p} Iv_t + m(C_0-C_i) + \frac{2Γ}{R} $$

$$ R = \frac {Γ/σ^*} {\frac{L}{c_p}P_t - 2 m P_c (1-k_0) C_i} $$

Given the following

$C_i = \frac{C_0}{1-(1-k_0)Iv_c}$ &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; - solute concentration of liquid at
 interface

$P_t = \frac{VR}{2α}$ &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; - thermal Péclet number

$P_c = \frac{VR}{2D}$ &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; - solutal Péclet number

$Iv_t(P_t) = P_t e^{P_t} E_1(P_t)$ &emsp;&emsp;&emsp;&emsp;- solutal Ivantsov function

$Iv_c(P_c) = P_c e^{P_c} E_1(P_c)$ &emsp;&emsp;&emsp;&ensp;&nbsp;- thermal Ivantsov function

$E_1(x) = \int_u^\infty \frac{exp(-s)}{s} ds$ &emsp;&emsp;&ensp;&emsp;&ensp;&nbsp;- first exponential
integral of $x$

The first equation calculates the LGK dendrite undercooling. It quantifies how much the bulk liquid must be cooled below
the equilibrium liquidus temperature at C0 to drive solidification. This extra undercooling allows 1. any accumulated temperature given out by the liquid solidifying to diffuse away from the dendrite tip, 2. allows reaching the lower
melting temperature caused by a build up of solute at the dendrite tip, and 3. allows overcoming the energy barrier
created by the surface energy of a high curvature dendrite tip. It uses dimensional analysis to solve for solute and
heat transport ahead of the parabaloid dendrite. Phase diagram constants are used to calculate the drop in liquidus temperature ahead of the solidification front due to solute enrichment.

The second equation calculates the LGK marginal stability criterion dendrite radius. A planar solidification front is
modified by adding a periodic pertubation. Too small and the curvature will drive the pertubation to shrink. Too large
and purtubation will grow by escaping the hot and solute rich solidification front. The dendrite radius is approximated
as the smallest pertubation that won't shrink. This gives an expression that is a function of the solute and temperature
field gradient ahead of the dendrite, which can be calculated for a parabaloid using the same dimensional analysis as
the first equation.

Note the extra factor of 2 in the second term of the second equation's denominator. Lipton, Glicksman, & Kurz remove
this factor in their paper in order to coerce this equation into agreeing with a prior published result for the case
where there is zero thermal field gradient and the second equation only depends on the solutal field gradient. This
change is not otherwise justified and is ignored in future iterations of this model such as LKT-BCT.

### The LKT-BCT Model

[J. Lipton, W. Kurz, R. Trivedi](https://doi.org/10.1016/0001-6160(87)90174-X) - [W.J. Boettinger, S.R. Coriell and R. 
Trivedi*](https://search.library.uq.edu.au/discovery/fulldisplay/alma991011497109703131/61UQ_INST:61UQ)

An extension of the LGK model that maintains accuracy at higher undercoolings and growth rates for alloys with linear
phase diagrams. 

$$ ∆T = \frac{L}{c_p} Iv_t + (mC_0 - m'C_i) + \frac{2Γ}{R} + \frac{V}{\mu} $$

$$ R = \frac {Γ/σ^*} {\frac{L}{c_p} P_t \xi_t  - 2 m P_c (1-k) C_i \xi_c} $$

Given the following:

$k = \frac{k_0 + (a_0V/D)}{1 + (a_0V/D)}$ &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; - velocity
 dependent partition coefficient

$C_i = \frac{C_0}{1-(1-k)Iv_c}$ &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; - solute
concentration of liquid at interface

$m' = m \left[1 + \frac{k_0 - k(1-ln(k/k_0))}{1-k_0} \right]$ &emsp;&emsp;&emsp;&emsp; - velocity dependent liquidus
slope

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

The derivation for $k$ assumes a small amount of solute, meaning the $-(1-k_0)\frac{C_0}{100}$ term that would otherwise
be in the denominator is negligible. While this assumption would be violated for higher $C_0$, especially when $k_0$ is
also small, including this term would mean that $k \not= k_0$ as $V \rightarrow 0$. As such, it is assumed that the
origional model authors were justified in keeping this assumption.

The first equation has an added fourth undercooling term, which specifies how much further the liquid must be cooled to
overcome the kinetic energy barrier for adding extra liquid atoms onto the solid lattice. This is assumed to be
negligible at the lower velocities expected in the LGK model. 

The LGK second equation assumes small Peclet numbers, where V*R << 1, meaning the stability functions ≈ 1. The second
LKT-BCT equation however drops this low V assumption, meaning these terms must be included.

The form of the equation above is exactly what BCT presented in their origional paper. More recently however, the
temperature used in $\mu$ for kinetic undercooling is often replaced with the bulk concentration dependent liquidus
temperature $ T_l = T_m + mC_0$, and the liquidus slope used in the second LKT-BCT equation is often replaced with the
non equilibrium velocity dependent liquidus slope $m'$. The latter form of this model arguably represents a solute rich,
fast moving interface more precisely.

\* While Boettinger, Coriell, and Trivedi's origional paper is not openly published online, it is famously well written
and the basis for this implementation. For a full derivation of this model, it can be accessed by requesting pages 13-25
of the linked conference paper. Note that a factor of $P_c$ is mistakenly missing in the denominator of the paper's R
equation.

### The CLW Model

[Chong-de Cao, Xiao-yu Lu, Bing-bo Wei](http://cpl.iphy.ac.cn/en/article/id/32758) 
 
An extension of LKT-BCT designed to better generalise to higher undercoolings and velocities for non-linear
phase diagrams, but makes strong assumptions and precise implementation details were never published. The below
implementation gives results close to those published.

$$ ∆T = \frac{L}{c_p} Iv_t + (mC_0 - m'C_i) + \frac{2Γ}{R} + \frac{V}{\mu} $$

$$ R = \frac {Γ/σ^*} {\frac{L}{c_p} P_t \xi_t  - 2 m P_c (1-k) C_i \xi_c} $$

Given the following:

$ T_0 = T_l(C_0) $ &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; - liquidus temperature at bulk
concentration

$ D = D_0 exp(-D_{Ea}/RT_{0}) $ &emsp;&emsp;&emsp;&emsp;&emsp; - temperature dependent diffusion coefficient of solute
in liquid

$ k0 = C_s(T_0) / C_l(T_0) $ &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; - temperature dependent equilibrium partition
coefficient

$ m = dT_l(C_i) / dC $ &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; - solute dependent liquidus slope

This model extends LKT-BCT to non-linear phase diagrams simply by using polynomial fits between the liquidus line,
solidus line, temperature, and solute concentration. While intuitive, many parts of the model still assume a linear
phase diagram. This includes the expression for the velocity dependant liquidus slope $m'$ used in the total
undercooling equation and bottom right term of the marginal stability criteria radius equation, where the liquidus
temperature gradient ahead of the dendrite tip is equal to the liquidus slope times the solute concentration gradient.

The papers based on this model also do not specify what temperatures each of the fit parameters should be evaluated at.
For a linear phase diagram $m$ and $k0$ are constant at all temperatures. For a non-linear phase diagram, the order at
which each undercooling is evaluated matters. Melt temperature changes between each undercooling component, so $m$, 
$k0$, $D$, and other derived parameters will also change. Here evaluating all parameters at $T_0$ gave the closest
result to those published, however the data points never match up exactly.
