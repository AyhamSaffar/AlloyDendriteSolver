# Alloys

Alloy structs are used to organise all thermodynamic constants for a given binary alloy system.

All *Alloy* objects can be used with the LGK model and must contain the following constants: 

- $L$ 	&emsp; Fusion enthalpy - $J/mol$
- $c_p$ &emsp; Melt heat capacity - $J/(mol K)$
- $m$ 	&emsp; Equilibrium liquidus slope - $K / C%$
- $k_0$	&emsp; Equilibrium partition coefficient - $C% / C%$
- $Γ$ 	&emsp; Gibbs-Thomson coefficient - $Km$
- $D$ 	&emsp; Diffusion coefficient of solute in liquid - $m^2/s$
- $α$	&emsp; Thermal diffusivity of liquid - $m^2/s$
- $σ^*$	&emsp; Stability constant - $m/m$

where C% is the concentration unit. Alloy varaible names with the *_wtp* suffix are in weight percent while the *_atp*
suffic refers to atom or mole percent.  

Optional parameters required for LKT-BCT models:

- $a_0$ &emsp; Solid atomic spacing - $m$
- $V_0$ &emsp; Speed of sound in liquid - $m/s$
- $T_m$ &emsp; Pure solid melting point at  - $K$

Further optional parameters required for CLW models, as well as the above:

- $D_0$ &emsp; Diffusivity at 0K - $m^2/s$
- $D_{Ea}$ &nbsp; Activation energy for diffusion - $J/mol$
- $T_l(C)$ &nbsp; Polynomial fit of concentration - $C%$ to liquidus temperature - $K$
- $C_s(T)$ &nbsp; Polynomial fit of temperature - $K$ to solidus concentration - $C%$ 
- $C_l(T)$ &nbsp; Polynomial fit of temperature - $K$ to liquidus concentration - $C%$

### Constants Varying With Other Parameters

Many of these constants actually vary with parameters such as temperature, solidification velocity, and solute mole
fraction. A model may incorporate some of these dependancies, but sometimes this can be safely ignored. These constants
may only vary a small amount in the range of parameters the model is designed for. Alternatively, the error created by
not including this variance may be negligible compared to the error created by other assumptions of a model.

### Choice of Units

The thermodynamic constants $L$ and $c_p$ have been normalised by mole throughout for consistency. Volumetric and
gravimetric values must be converted to molar for use in this library.

Concentration based quantities ($m$ and $k_0$) can either be in at.% or wt.%. Care must be taken to ensure all
concentration values have the same units as the corresponding Alloy.

Units for $k_0$ are often omitted as it is unitless overall. However it has a different value when in wt.% / wt.%
compared to at.% / at.%. Consider the situation where the solute metal atoms have a much lower atomic mass then
the bulk metal. If The solute has three times the wt.% in the liquid as the solid ($k_0$ = 0.33 wt.%/wt.%), it follows
that it would have more than three times the atom fraction in the liquid as the solid. These units are easy to miss and
cannot be easily interconverted in most cases.
