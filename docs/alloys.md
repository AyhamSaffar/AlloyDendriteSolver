# Alloys

Alloy structs are used to organise all thermodynamic constants for a given binary alloy system.

All *Alloy* objects must contain the following constants: 

- $L$ 	&nbsp; Fusion enthalpy - $J/mol$
- $c_p$ &nbsp; Melt heat capacity - $J/(mol K)$
- $m$ 	&nbsp; Equilibrium liquidus slope - $K / C\%$
- $k_0$	&nbsp; Equilibrium partition coefficient - $C\% / C\%$
- $Γ$ 	&nbsp; Gibbs-Thomson coefficient - $Km$
- $D$ 	&nbsp; Diffusion coefficient of solute in liquid - $m^2/s$
- $α$	&nbsp; Thermal diffusivity of liquid - $m^2/s$
- $σ^*$	&nbsp; Stability constant - $m/m$

where C% is the concentration unit. Alloy varaible names with the *_wtp* suffix are in weight percent while the *_atp*
suffic refers to atom or mole percent.  

Optional parameters required for LKT-BCT models:

- $a_0$ &nbsp; Solid atomic spacing - $m$
- $V_0$ &nbsp; Speed of sound in liquid - $m/s$
- $T_m$ &nbsp; Pure solid melting point at  - $K$

### Constants Varying With Other Parameters

Many of these constants actually vary with parameters such as temperature, solidification velocity, and solute mole
fraction. A model may incorporate some of these dependancies, but sometimes this can be safely ignored. These constants
may only vary a small amount in the range of parameters the model is designed for. Alternatively, the error created by
not including this variance may be negligible compared to the error created by other assumptions of a model.

### Choice of Units

The thermodynamic constants $L$ and $c_p$ have been normalised by mole throughout for consistency. Volumetric and
gravimetric values must be converted to molar for use in this library.

Concentration based quantities ($m$ and $k_0$) can either be in at.% or wt.%. Alloys only used with the LGK model can be
either as this model only uses these quantities in their unitless forms during dimensional analysis calculations.
**Alloys used with the LKT-BCT model must be in at.%** as this model contains terms that require the proportion of atoms
in a given phase.

Units for $k_0$ are often omitted as it is unitless overall. However it has a different value when in wt.% / wt.%
compared to at.% / at.%. Consider the situation where the solute metal atoms have a much lower atomic mass then
the bulk metal. If The solute has three times the wt.% in the liquid as the solid ($k_0$ = 0.33 wt.%/wt.%), it follows
that it would have more than three times the atom fraction in the liquid as the solid. These units are easy to miss and
cannot be easily interconverted in most cases. A full phase diagram line fit is required.
