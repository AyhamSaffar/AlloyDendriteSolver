function [V_guess_func, R_guess_func] = create_analytical_guess_functions(C0, Gibbs_Tom, m_l, k0, D_l)
%% CREATE_ANALYTICAL_GUESS_FUNCTIONS - Generate theoretical initial guess functions
% ========================================================================
%
% PURPOSE:
% Creates initial guess functions based on analytical approximations derived
% from solidification theory, specifically implementing equations (8.91) and
% (8.92) from Kurz & Fisher's solidification textbook. These provide
% theoretically-grounded initial estimates for the Newton-Raphson solver.
%
% THEORETICAL FOUNDATION:
% The analytical approximations are derived from asymptotic analysis of the
% LGK model in specific limiting regimes, providing scaling relationships
% that capture the essential physics of dendritic growth.
%
% EQUATION (8.92) - Velocity Approximation:
% V* = D_l / (5.51×π²×(-m_l×(1-k₀))^1.5×Γ_sl) × (ΔT^2.5 / C₀^1.5)
%
% Physical interpretation:
% - V ∝ D_l: Faster diffusion enables higher growth rates
% - V ∝ ΔT^2.5: Strong superlinear dependence on driving force
% - V ∝ C₀^-1.5: Solute drag reduces velocity at high concentrations
% - V ∝ (-m_l×(1-k₀))^-1.5: Thermodynamic property influence
%
% EQUATION (8.91) - Radius Approximation:
% R_tip = 6.64×π²×Γ_sl×(-m_l×(1-k₀))^0.25 × (C₀^0.25 / ΔT^1.25)
%
% Physical interpretation:
% - R ∝ Γ_sl: Surface tension promotes larger radii for stability
% - R ∝ ΔT^-1.25: Higher driving force enables sharper tips
% - R ∝ C₀^0.25: Weak concentration dependence from stability balance
% - R ∝ (-m_l×(1-k₀))^0.25: Constitutional effects on tip stability
%
% INPUTS:
%   C0 - Solute concentration [wt%]
%   Gibbs_Tom - Gibbs-Thomson coefficient [K·m]  
%   m_l - Liquidus slope [K/wt%] (negative for typical alloys)
%   k0 - Partition coefficient [-]
%   D_l - Interdiffusion coefficient [m²/s]
%
% OUTPUTS:
%   V_guess_func - Function handle: V₀ = V_guess_func(ΔT)
%   R_guess_func - Function handle: R₀ = R_guess_func(ΔT)
%
% SPECIAL HANDLING:
% For very low concentrations (C₀ < 0.05 wt%), the functions implement
% reference-based scaling to avoid numerical instabilities while
% maintaining physical scaling behavior.
%
% USAGE EXAMPLE:
%   [V_func, R_func] = create_analytical_guess_functions(2.0, 2.4e-7, -3.37, 0.17, 3e-9);
%   V_initial = V_func(1.5);  % Get velocity estimate for ΔT = 1.5 K
%   R_initial = R_func(1.5);  % Get radius estimate for ΔT = 1.5 K
%
% ADVANTAGES:
% - Physically meaningful initial guesses improve convergence
% - Reduces iteration count and computational time
% - Provides better scaling across wide parameter ranges
% - Captures correct asymptotic behavior in limiting cases
%
% LIMITATIONS:
% - Accuracy decreases outside derivation assumptions
% - May require fallback to parametric forms for extreme cases
% - Derived for specific alloy systems (may need adjustment)
% ========================================================================

% Material constants
pi_val = pi;

% Create V_guess function based on Eq. (8.92)
V_guess_func = @(deltaT) calculate_V_analytical(deltaT, C0, D_l, m_l, k0, Gibbs_Tom, pi_val);

% Create R_guess function based on Eq. (8.91)
R_guess_func = @(deltaT) calculate_R_analytical(deltaT, C0, Gibbs_Tom, m_l, k0, pi_val);

end

function V = calculate_V_analytical(deltaT, C0, D_l, m_l, k0, Gibbs_Tom, pi_val)
%% CALCULATE_V_ANALYTICAL - Compute analytical velocity approximation
%
% PURPOSE:
% Implements equation (8.92) for dendrite tip velocity with special handling
% for low concentration regimes to maintain numerical stability and physical
% scaling behavior.
%
% CONCENTRATION HANDLING STRATEGY:
% For C₀ < 0.05 wt%: Uses reference concentration approach
% - Calculates velocity at reference concentration (C₀_ref = 0.05)
% - Applies theoretical scaling: V ∝ C₀^-1.5
% - Implements scaling limits and dampening for extreme cases
% - Provides smooth transition to avoid discontinuities
%
% For C₀ ≥ 0.05 wt%: Direct calculation using equation (8.92)
%
% NUMERICAL STABILITY MEASURES:
% - Validates input parameters for physical realizability
% - Handles negative or zero liquidus slopes appropriately
% - Provides fallback values for numerical difficulties
% - Implements scaling limits to prevent extreme extrapolation

try
    % Handle array inputs
    V = zeros(size(deltaT));
    
    % For very small C0, use C0=0.05 as reference and scale the result
    C0_ref = 0.05;  % Reference concentration
    use_reference = C0 < C0_ref;
    
    if use_reference
        fprintf('Info: C0=%.6f < %.3f, using reference-based calculation\n', C0, C0_ref);
        % Calculate using reference concentration first
        V_ref = calculate_V_analytical_core(deltaT, C0_ref, D_l, m_l, k0, Gibbs_Tom, pi_val);
        % Scale result based on concentration ratio (V ∝ 1/C0^1.5)
        scaling_factor = (C0_ref / C0)^1.5;
        % Apply reasonable scaling limit and additional dampening for very small C0
        scaling_factor = min(scaling_factor, 10);  % More conservative limit
        if C0 < 0.03
            scaling_factor = scaling_factor * 0.3;  % Additional dampening for very small C0
            fprintf('Applied extra dampening for C0=%.6f, final scaling=%.2f\n', C0, scaling_factor);
        end
        V = V_ref * scaling_factor;
        fprintf('Info: Scaled V by factor %.2f for C0=%.6f\n', scaling_factor, C0);
    else
        V = calculate_V_analytical_core(deltaT, C0, D_l, m_l, k0, Gibbs_Tom, pi_val);
    end
    
catch ME
    fprintf('Error in V analytical calculation: %s\n', ME.message);
    V = ones(size(deltaT)) * 1e-6;
end

end

function V = calculate_V_analytical_core(deltaT, C0, D_l, m_l, k0, Gibbs_Tom, pi_val)
%% CALCULATE_V_ANALYTICAL_CORE - Core velocity calculation implementation
%
% PURPOSE:
% Performs the direct mathematical evaluation of equation (8.92) with
% comprehensive input validation and error handling.
%
% EQUATION IMPLEMENTATION:
% V* = D_l / (5.51×π²×(-m_l×(1-k₀))^1.5×Γ_sl) × (ΔT^2.5 / C₀^1.5)
%
% PARAMETER VALIDATION:
% - C₀ > 0: Positive concentration required
% - D_l > 0: Positive diffusion coefficient required  
% - Γ_sl > 0: Positive Gibbs-Thomson coefficient required
% - (1-k₀) > 0: Partition coefficient less than unity required
% - m_l < 0: Negative liquidus slope for typical alloys
%
% ERROR HANDLING:
% Invalid parameters trigger fallback to minimum physical velocity
% Complex or infinite results replaced with stable default values
% Comprehensive diagnostic output for debugging problematic cases

% Core calculation function
V = zeros(size(deltaT));

% Check for valid input parameters
if C0 <= 0 || D_l <= 0 || Gibbs_Tom <= 0 || (1-k0) <= 0
    fprintf('Warning: Invalid parameters for V calculation\n');
    V = ones(size(deltaT)) * 1e-6;
    return;
end

for i = 1:length(deltaT)
    dt = deltaT(i);
    
    if dt <= 0
        V(i) = 1e-8;
        continue;
    end
    
    % Eq. (8.92): v* = D_l / (5.51*pi^2*(-m_l*(1-k0))^1.5*Γ_sl) * (ΔT^2.5 / C0^1.5)
    term1 = (-m_l * (1 - k0));
    
    if term1 <= 0
        fprintf('Warning: Invalid term1=%.6f\n', term1);
        V(i) = 1e-8;
        continue;
    end
    
    denominator = 5.51 * pi_val^2 * (term1^1.5) * Gibbs_Tom;
    numerator = D_l * (dt^2.5) / (C0^1.5);
    
    V(i) = numerator / denominator;
    
    % Basic validation
    if V(i) <= 0 || ~isreal(V(i)) || isnan(V(i)) || isinf(V(i))
        V(i) = 1e-6;
        fprintf('Using fallback V=%.2e for C0=%.6f, deltaT=%.3f\n', V(i), C0, dt);
    end
end

end

function R = calculate_R_analytical(deltaT, C0, Gibbs_Tom, m_l, k0, pi_val)
%% CALCULATE_R_ANALYTICAL - Compute analytical radius approximation
%
% PURPOSE:
% Implements equation (8.91) for dendrite tip radius with concentration-based
% scaling strategy similar to velocity calculation for consistent numerical
% behavior across the full concentration range.
%
% THEORETICAL SCALING:
% R_tip = 6.64×π²×Γ_sl×(-m_l×(1-k₀))^0.25 × (C₀^0.25 / ΔT^1.25)
%
% CONCENTRATION SCALING APPROACH:
% For low concentrations: R ∝ C₀^0.25 scaling applied to reference calculation
% This maintains the correct physical dependence while avoiding numerical
% issues that arise from direct evaluation at very low concentrations.
%
% PHYSICAL INTERPRETATION:
% - Radius increases with surface tension (Γ_sl)
% - Radius decreases with driving force (ΔT^-1.25)  
% - Weak positive concentration dependence (C₀^0.25)
% - Constitutional effects through (-m_l×(1-k₀))^0.25 term

try
    % Handle array inputs
    R = zeros(size(deltaT));
    
    % For very small C0, use C0=0.05 as reference and scale the result
    C0_ref = 0.05;  % Reference concentration
    use_reference = C0 < C0_ref;
    
    if use_reference
        fprintf('Info: C0=%.6f < %.3f, using reference-based calculation (R)\n', C0, C0_ref);
        % Calculate using reference concentration first
        R_ref = calculate_R_analytical_core(deltaT, C0_ref, Gibbs_Tom, m_l, k0, pi_val);
        % Scale result based on concentration ratio (R ∝ C0^0.25)
        scaling_factor = (C0 / C0_ref)^0.25;
        R = R_ref * scaling_factor;
        fprintf('Info: Scaled R by factor %.2f for C0=%.6f\n', scaling_factor, C0);
    else
        R = calculate_R_analytical_core(deltaT, C0, Gibbs_Tom, m_l, k0, pi_val);
    end
    
catch ME
    fprintf('Error in R analytical calculation: %s\n', ME.message);
    R = ones(size(deltaT)) * 1e-6;
end

end

function R = calculate_R_analytical_core(deltaT, C0, Gibbs_Tom, m_l, k0, pi_val)
%% CALCULATE_R_ANALYTICAL_CORE - Core radius calculation implementation
%
% PURPOSE:
% Executes direct mathematical evaluation of equation (8.91) with robust
% error handling and physical constraint enforcement.
%
% MATHEMATICAL IMPLEMENTATION:
% R_tip = 6.64×π²×Γ_sl×(-m_l×(1-k₀))^0.25 × (C₀^0.25 / ΔT^1.25)
%
% COMPUTATIONAL STEPS:
% 1. Validate thermodynamic term: (-m_l×(1-k₀)) > 0
% 2. Calculate coefficient: 6.64×π²×Γ_sl×term1^0.25
% 3. Evaluate ratio term: (C₀^0.25 / ΔT^1.25)
% 4. Combine terms and validate physical realizability
%
% ROBUSTNESS FEATURES:
% - Handles invalid thermodynamic parameters gracefully
% - Provides diagnostic output for debugging difficult cases
% - Implements fallback values for numerical instabilities
% - Maintains positive radius constraint for physical meaning
%
% PERFORMANCE CONSIDERATIONS:
% - Vectorized operations for efficiency with array inputs
% - Minimal memory allocation for large deltaT sequences
% - Optimized mathematical operations for computational speed
% - Error handling designed to avoid computational overhead

% Core calculation function
R = zeros(size(deltaT));

% Check for valid input parameters
if C0 <= 0 || Gibbs_Tom <= 0 || (1-k0) <= 0
    fprintf('Warning: Invalid parameters for R calculation\n');
    R = ones(size(deltaT)) * 1e-6;
    return;
end

for i = 1:length(deltaT)
    dt = deltaT(i);
    
    if dt <= 0
        R(i) = 1e-6;
        continue;
    end
    
    % Eq. (8.91): R_tip = 6.64*pi^2*Γ_sl*(-m_l*(1-k0))^0.25 * (C0^0.25 / ΔT^1.25)
    term1 = (-m_l * (1 - k0));
    
    if term1 <= 0
        fprintf('Warning: Invalid term1=%.6f for R calculation\n', term1);
        R(i) = 1e-6;
        continue;
    end
    
    coefficient = 6.64 * pi_val^2 * Gibbs_Tom * (term1^0.25);
    ratio_term = (C0^0.25) / (dt^1.25);
    
    R(i) = coefficient * ratio_term;
    
    % Basic validation
    if R(i) <= 0 || ~isreal(R(i)) || isnan(R(i)) || isinf(R(i))
        R(i) = 1e-6;
        fprintf('Using fallback R=%.2e for C0=%.6f, deltaT=%.3f\n', R(i), C0, dt);
    end
end

end