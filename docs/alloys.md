# Alloys

Alloy structs are used to organise all thermodynamic constants for a given binary alloy system.

All *Alloy* object must contain the following constants: 

- $L$ 	&nbsp; Fusion enthalpy - $J/mol$
- $c_p$ &nbsp; Melt heat capacity - $J/(mol K)$
- $m$ 	&nbsp; Equilibrium liquidus slope - $K / wt.\%$
- $k_0$	&nbsp; Partition coefficient - $wt.\% / wt.\%$
- $Γ$ 	&nbsp; Gibbs-Thomson coefficient - $Km$
- $D$ 	&nbsp; Diffusion coefficient of solute in liquid - $m^2/s$
- $α$	&nbsp; Thermal conductivity of liquid - $m^2/s$
- $σ^*$	&nbsp; Stability constant - $m/m$

Optional parameters required for some models:

- $a_0$ &nbsp; Solid atomic spacing - $m$
- $V_0$ &nbsp; Speed of sound in liquid - $m/s$
- $T_m$ &nbsp; Solid melting point - $K$

### Constants Varying With Other Parameters

Many of these constants actually vary with parameters such as temperature, solidification velocity, and solute mole
fraction. A model may incorporate some of these dependancies, but sometimes this can be safely ignored. These constants
may only vary a small amount in the range of parameters the model is designed for. Alternatively, the error created by
not including this variance may be negligible compared to the error created by other assumptions of a model.

### Choice of Units

Thermodynamic constants have been normalised by mole while phase diagram constants have been normalised by weight
percent. The choice of units often do not matter as these constants appear as ratios (E.G $L/C_p$) or are converted to
unitless forms (E.G. during dimensional analysis calculations). However some models use thermodynamic constans on their
own (E.G. $L$ in the kinetic undercooling term in LKT_BCT), so molar normalisation avoids the need for easy to miss
atomic mass terms.

Units for $k_0$ are often omitted as it is unitless overall. However it has a different value when in wt.% / wt.%
compared to mol frac / mol frac. Consider the situation where the solute metal atoms have a much lower atomic mass then
the bulk metal. If The solute has three times the wt.% in the liquid as the solid ($k_0$ = 0.33 wt.%/wt.%), it follows
that it would have more than three times the atom fraction in the liquid as the solid. This conversion is easy to miss
and so $k_0$ units should be checked.
