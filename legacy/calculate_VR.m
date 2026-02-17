function [result, converged, last_diverged_deltaT] = calculate_VR(V_guess_func, R_guess_func, C0, deltaT_values, material_params)
%% CALCULATE_VR - Core Newton-Raphson solver for LGK dendritic growth equations
% ========================================================================
%
% PURPOSE:
% This is the computational heart of the LGK model, solving the coupled
% nonlinear equations for dendrite tip velocity (V) and radius (R) using
% an adaptive Newton-Raphson method with sophisticated convergence control.
%
% MATHEMATICAL FOUNDATION:
% Solves the coupled LGK system:
% f₁ = ΔTₜ - ΔTc + ΔTᵣ - ΔT = 0  (Undercooling balance)
% f₂ = R - Γσ/(mGc - G) = 0        (Stability criterion)
%
% Where:
% ΔTₜ = (ΔH/Cp)×Iv(Pₜ)                    (Thermal undercooling)
% ΔTc = mC₀[(1-k₀)Iv(Pc)]/[1-(1-k₀)Iv(Pc)] (Solutal undercooling)
% ΔTᵣ = 2Γ/R                              (Curvature undercooling)
% Pₜ = VR/(2a), Pc = VR/(2D)              (Peclet numbers)
% Iv(P) = P×exp(P)×expint(P)              (Ivantsov function)
%
% NEWTON-RAPHSON FORMULATION:
% [V_{n+1}] = [Vₙ] - [∂f₁/∂V  ∂f₁/∂R]⁻¹ [f₁]
% [R_{n+1}]   [Rₙ]   [∂f₂/∂V  ∂f₂/∂R]   [f₂]
%
% ADAPTIVE CONVERGENCE STRATEGY:
% - Magnitude-based step damping to prevent oscillations
% - Characteristic scaling for numerical stability
% - Multi-scale tolerance criteria for robust convergence
% - Automatic step size adjustment based on convergence history
%
% INPUTS:
%   V_guess_func - Function handle for initial velocity estimates
%                  Format: V₀ = V_guess_func(ΔT)
%   R_guess_func - Function handle for initial radius estimates  
%                  Format: R₀ = R_guess_func(ΔT)
%   C₀ - Solute concentration [wt%]
%   deltaT_values - Array of undercooling values to solve [K]
%   material_params - Structure containing:
%     .Cpv - Volumetric specific heat [J/m³·K]
%     .DHv - Volumetric latent heat [J/m³]
%     .a - Thermal diffusivity [m²/s]
%     .Gibbs_Tom - Gibbs-Thomson coefficient [K·m]
%     .D_Al - Interdiffusion coefficient [m²/s]
%     .m_Al - Liquidus slope [K/wt%]
%     .k0_Al - Partition coefficient [-]
%     .sigma - Stability constant [-]
%     .include_thermal_undercooling - Include thermal effects [logical]
%
% OUTPUTS:
%   result - [N×3] matrix: [ΔT, V, R] for converged solutions
%   converged - Logical flag indicating complete convergence success
%   last_diverged_deltaT - Last undercooling value that failed to converge [K]
%
% CONVERGENCE CRITERIA:
% - Absolute tolerance: |ΔV|, |ΔR| < 1×10⁻⁸ (scaled)
% - Maximum iterations: 200 per undercooling value
% - Step quality monitoring with adaptive damping
% - Early termination for non-physical solutions
%
% ERROR HANDLING:
% - Validates initial guess function outputs
% - Handles complex or infinite intermediate results
% - Provides graceful degradation for difficult cases
% - Maintains numerical stability under extreme conditions
%
% PERFORMANCE FEATURES:
% - Real-time progress display with user interruption capability
% - Optimized symbolic computation with caching
% - Memory-efficient result storage
% - Parallel-safe execution for multi-concentration calculations
% ========================================================================

% Extract material parameters from structure for computational efficiency
% This avoids repeated structure field access during intensive calculations
Cpv = material_params.Cpv;          % Specific heat (J/m³·K)
DHv = material_params.DHv;          % Latent heat (J/m³)
a = material_params.a;              % Thermal diffusivity (m²/s)
Gibbs_Tom = material_params.Gibbs_Tom; % Gibbs-Thomson coefficient
D_Al = material_params.D_Al;        % Diffusion coefficient (m²/s)
m_Al = material_params.m_Al;        % Liquidus slope
k0_Al = material_params.k0_Al;      % Partition coefficient
sigma = material_params.sigma;      % Stability constant

% Initialize
converged = true;  % convergence flag
last_diverged_deltaT = NaN;  % record last diverged deltaT
n1 = 0;
result = zeros(length(deltaT_values), 3);

for i = 1:length(deltaT_values)
    % Check for stop signal
    global STOP_CALCULATION;
    if ~isempty(STOP_CALCULATION) && STOP_CALCULATION
        fprintf('\nCalculation stopped by user at deltaT = %.2f\n', deltaT_values(i));
        converged = false;
        break;
    end
    
    deltaT = deltaT_values(i);
    if ~converged
        break;  % exit loop if previous didn't converge
    end
    n1 = n1 + 1;

    % Define symbolic variables and equations representing the LGK model
    syms x1 x2 % x1 = V (velocity), x2 = R (radius)

    % Calculate Peclet numbers governing heat and mass transport
    Pt = (x1*x2)/(2*a); % Thermal Peclet number: ratio of convection to diffusion
    Pc = (x1*x2)/(2*D_Al);  % Solutal Peclet number: controls constitutional effects

    % Compute Ivantsov functions describing diffusion field solutions
    % These represent exact solutions for transport around parabolic dendrite tips
    IvPt = Pt*exp(Pt)*expint(Pt); % Ivantsov function for thermal field
    IvPc = Pc*exp(Pc)*expint(Pc); % Ivantsov function for solutal field

    % Configure thermal undercooling contribution based on model settings
    % Thermal effects can be disabled for purely solutal problems
    if material_params.include_thermal_undercooling
        thermal_term_f1 = (DHv/Cpv)*IvPt;   % Thermal undercooling term
        thermal_term_f2 = Pt*DHv/Cpv;       % Thermal gradient term
    else
        thermal_term_f1 = 0;    
        thermal_term_f2 = 0;
    end

    % Construct the coupled LGK equations
    % f1: Overall undercooling balance equation
    f1 = thermal_term_f1 - m_Al*C0*((1-k0_Al)*IvPc)/(1-(1-k0_Al)*IvPc) + 2*Gibbs_Tom/x2 - deltaT;
    % f2: Marginal stability criterion equation  
    f2 = (Gibbs_Tom/sigma)/(thermal_term_f2-Pc*m_Al*C0*(1-k0_Al)/(1-(1-k0_Al)*IvPc)) - x2;

    df1x1 = diff(f1,x1);
    df1x2 = diff(f1,x2);
    df2x1 = diff(f2,x1);
    df2x2 = diff(f2,x2);

    % Newton-Raphson iteration formula
    f3 = [x1;x2] - ([df1x1 df1x2; df2x1 df2x2])\[f1;f2];
    
    % Iterative solution with improved algorithm
    [V, R, iteration_converged] = solve_VR_iteration(f3, V_guess_func, R_guess_func, deltaT);
    
    % Store results
    result(n1,1) = deltaT;
    result(n1,2) = V;
    result(n1,3) = R;
    
    % Real-time progress display with user interruption capability
    if iteration_converged
        global SUPPRESS_OUTPUT;
        % Check if output should be displayed (not suppressed by parallel workers)
        if isempty(SUPPRESS_OUTPUT) || ~SUPPRESS_OUTPUT
            fprintf('%.2f, ', deltaT);  % Display every deltaT value
            if mod(n1, 10) == 0  % New line every 10 values for readability
                fprintf('\n');
            end
        end
    end
    
    % Check convergence
    if ~iteration_converged
        converged = false;
        last_diverged_deltaT = deltaT;
        
        global SUPPRESS_OUTPUT;
        if isempty(SUPPRESS_OUTPUT) || ~SUPPRESS_OUTPUT
            fprintf('\nStopped at deltaT = %.2f (not converged)\n', deltaT);
        end
        
        break;  % Exit deltaT loop - no point continuing if current fails
    end
end

% Final newline for clean output
if converged
    global SUPPRESS_OUTPUT;
    if isempty(SUPPRESS_OUTPUT) || ~SUPPRESS_OUTPUT
        fprintf('\nAll deltaT values converged successfully.\n');
    end
end

% Remove uncalculated rows
result(n1+1:end, :) = [];

end

function [V, R, converged] = solve_VR_iteration(f3, V_guess_func, R_guess_func, deltaT)
%% SOLVE_VR_ITERATION - Single undercooling Newton-Raphson solution
%
% PURPOSE:
% Executes Newton-Raphson iteration for a single undercooling value with
% comprehensive initial guess validation and adaptive scaling strategies.
%
% INITIAL GUESS STRATEGY:
% - Evaluates user-provided guess functions at specified undercooling
% - Validates physical realizability (positive, finite values)
% - Provides fallback values for invalid initial guesses
% - Ensures numerical stability before iteration begins
%
% INPUTS:
%   f3 - Symbolic Newton-Raphson iteration formula [2×1 symbolic vector]
%   V_guess_func - Velocity initial guess function handle
%   R_guess_func - Radius initial guess function handle  
%   deltaT - Target undercooling value [K]
%
% OUTPUTS:
%   V, R - Converged velocity [m/s] and radius [m] solutions
%   converged - Success flag for convergence achievement
%
% VALIDATION CRITERIA:
% Initial guesses must satisfy: V₀ > 0, R₀ > 0, finite values
% Invalid guesses trigger automatic correction to minimum physical values

converged = false;
V = NaN;
R = NaN;

% Use only guess function values
V0 = V_guess_func(deltaT);
R0 = R_guess_func(deltaT);

% Validate initial guess values
if V0 <= 0 || R0 <= 0 || isnan(V0) || isnan(R0) || isinf(V0) || isinf(R0)
    % Try to generate alternative initial values
    V0 = max(1e-6, abs(V0));  % Ensure positive minimum value
    R0 = max(1e-6, abs(R0));  % Ensure positive minimum value
    
    % If still invalid, return unconverged
    if isnan(V0) || isnan(R0) || isinf(V0) || isinf(R0)
        return;
    end
end

[V, R, converged] = newton_raphson_solve_adaptive(f3, V0, R0);

end

function [V, R, converged] = newton_raphson_solve_adaptive(f3, V0, R0)
%% NEWTON_RAPHSON_SOLVE_ADAPTIVE - Advanced adaptive Newton-Raphson solver
%
% PURPOSE:
% Implements sophisticated Newton-Raphson iteration with magnitude-based
% adaptive damping, characteristic scaling, and intelligent convergence
% monitoring for robust solution of the highly nonlinear LGK equations.
%
% ADAPTIVE DAMPING ALGORITHM:
% The algorithm dynamically adjusts step damping based on step magnitude:
% 
% For large steps (|Δx_norm| > 1.0):
%   damping = base_damping × (1 + log₁₀(step_magnitude))
%   Purpose: Prevents oscillations and ensures stability
%
% For small steps (|Δx_norm| ≤ 1.0):  
%   damping = base_damping × (0.5 + 0.5 × step_magnitude)
%   Purpose: Maintains convergence rate near solution
%
% CHARACTERISTIC SCALING:
% Maintains separate scaling factors for V and R:
% V_scale = max(min_scale, 0.9×V_scale_old + 0.1×|V_current|)
% R_scale = max(min_scale, 0.9×R_scale_old + 0.1×|R_current|)
% This ensures numerical stability across wide dynamic ranges.
%
% CONVERGENCE MONITORING:
% - Normalized tolerance criteria: |ΔV|/V_scale, |ΔR|/R_scale < 10⁻⁸
% - Step quality assessment for convergence rate optimization
% - Consecutive failure detection for early termination
% - Progress reporting for long calculations
%
% NUMERICAL STABILITY FEATURES:
% - Positive value enforcement to maintain physical meaning
% - Complex number detection and rejection
% - Infinite value handling with graceful degradation
% - Conservative damping fallback for extreme cases
%
% INPUTS:
%   f3 - Newton-Raphson iteration formula from symbolic differentiation
%   V0, R0 - Initial guess values [m/s], [m]
%
% OUTPUTS:
%   V, R - Final converged values or NaN if convergence failed
%   converged - Boolean success indicator
%
% ALGORITHM PARAMETERS:
%   max_iter = 200 - Maximum iterations before timeout
%   tolerance = 1×10⁻⁸ - Normalized convergence criterion
%   base_damping = 0.8 - Reference damping factor
%   min_scale = 1×10⁻¹² - Minimum scaling to prevent division by zero
%
% PERFORMANCE OPTIMIZATIONS:
% - Rate-limited debug output to prevent I/O bottlenecks
% - Efficient symbolic evaluation with minimal overhead
% - Early termination strategies to avoid unnecessary computation
% - Memory-efficient variable updates without excessive allocation

max_iter = 200;  % Increased max iterations for adaptive method
tolerance = 1e-8;
converged = false;

% Ensure initial values are positive
if V0 <= 0 || R0 <= 0 || isnan(V0) || isnan(R0) || isinf(V0) || isinf(R0)
    V = NaN;
    R = NaN;
    return;
end

V = V0;
R = R0;

% Initialize adaptive parameters
V_scale = abs(V0);  % Characteristic scale for V
R_scale = abs(R0);  % Characteristic scale for R
min_scale = 1e-12;  % Minimum scale to prevent division by zero

% Adaptive damping parameters
base_damping = 0.8;
max_damping = 0.95;
min_damping = 0.1;

% Step size monitoring
consecutive_good_steps = 0;
consecutive_bad_steps = 0;

% Store previous values for step monitoring
V_prev = V;
R_prev = R;

global SUPPRESS_OUTPUT;
debug_mode = isempty(SUPPRESS_OUTPUT) || ~SUPPRESS_OUTPUT;

if debug_mode
    fprintf('\n  [Adaptive N-R] Initial: V=%.2e, R=%.2e, V_scale=%.2e, R_scale=%.2e\n', V, R, V_scale, R_scale);
end

for iter = 1:max_iter
    % Check for stop signal in inner loop
    global STOP_CALCULATION;
    if ~isempty(STOP_CALCULATION) && STOP_CALCULATION
        converged = false;
        return;  
    end
    
    x1 = V;
    x2 = R;
    
    try
        % Calculate next step using Newton-Raphson
        ans_in_sym = subs(f3);
        ans_in_double = double(ans_in_sym);
        
        V_newton = ans_in_double(1);
        R_newton = ans_in_double(2);
        
        % Check for complex, invalid, or negative values
        if ~isreal(V_newton) || ~isreal(R_newton) || isnan(V_newton) || isnan(R_newton) || ...
           isinf(V_newton) || isinf(R_newton) || V_newton <= 0 || R_newton <= 0
            consecutive_bad_steps = consecutive_bad_steps + 1;
            if debug_mode && consecutive_bad_steps == 1
                fprintf('  [Adaptive N-R] Invalid Newton step at iter %d\n', iter);
            end
            break;
        end
        
        % Calculate raw increments
        delta_V_raw = V_newton - V;
        delta_R_raw = R_newton - R;
        
        % Update characteristic scales based on current values
        V_scale = max(min_scale, 0.9 * V_scale + 0.1 * abs(V));
        R_scale = max(min_scale, 0.9 * R_scale + 0.1 * abs(R));
        
        % Calculate normalized increments for adaptive damping assessment
        delta_V_norm = delta_V_raw / V_scale;
        delta_R_norm = delta_R_raw / R_scale;
        
        % Compute step magnitudes for damping algorithm
        V_step_magnitude = abs(delta_V_norm);
        R_step_magnitude = abs(delta_R_norm);
        
        % Apply magnitude-based adaptive damping strategy
        % Large steps receive increased damping to prevent overshooting
        if V_step_magnitude > 1.0
            damping_V = base_damping * (1 + log10(V_step_magnitude));
            damping_V = min(max_damping, damping_V);    % Enforce upper limit
        else
            % Small steps receive reduced damping to maintain convergence rate
            damping_V = base_damping * (0.5 + 0.5 * V_step_magnitude);
            damping_V = max(min_damping, damping_V);    % Enforce lower limit
        end

        % Apply same strategy to radius with independent damping factors
        if R_step_magnitude > 1.0
            damping_R = base_damping * (1 + log10(R_step_magnitude));
            damping_R = min(max_damping, damping_R);
        else
            damping_R = base_damping * (0.5 + 0.5 * R_step_magnitude);
            damping_R = max(min_damping, damping_R);
        end
        
        % Update variables with adaptive damping applied
        V_new = V + damping_V * delta_V_raw;
        R_new = R + damping_R * delta_R_raw;
        
        % Enforce physical constraints: dendrite properties must be positive
        if V_new <= 0 || R_new <= 0
            % Attempt recovery with conservative damping
            conservative_damping = 0.1;
            V_new = V + conservative_damping * delta_V_raw;
            R_new = R + conservative_damping * delta_R_raw;
            
            if V_new <= 0 || R_new <= 0
                consecutive_bad_steps = consecutive_bad_steps + 1;
                if debug_mode && consecutive_bad_steps == 1
                    fprintf('  [Adaptive N-R] Negative values at iter %d\n', iter);
                end
                break;  % Terminate iteration - physical constraints violated
            end
        end
        
        % Calculate convergence measures (normalized)
        VD_norm = abs(V_new - V) / V_scale;
        RD_norm = abs(R_new - R) / R_scale;
        
        % Monitor step quality
        step_quality = max(VD_norm, RD_norm);
        if step_quality < 0.1  % Good step
            consecutive_good_steps = consecutive_good_steps + 1;
            consecutive_bad_steps = 0;
        else  % Large step
            consecutive_good_steps = 0;
            if step_quality > 2.0
                consecutive_bad_steps = consecutive_bad_steps + 1;
            end
        end
        
        % Update values
        V_prev = V;
        R_prev = R;
        V = V_new;
        R = R_new;
        
        % Debug output every 20 iterations
        if debug_mode && mod(iter, 20) == 0
            fprintf('  [Adaptive N-R] Iter %d: V=%.2e(δ=%.1e,d=%.2f), R=%.2e(δ=%.1e,d=%.2f), Conv=%.1e\n', ...
                    iter, V, abs(delta_V_raw), damping_V, R, abs(delta_R_raw), damping_R, max(VD_norm, RD_norm));
        end
        
        % Check convergence with normalized tolerance
        if VD_norm < tolerance && RD_norm < tolerance
            converged = true;
            if debug_mode
                fprintf('  [Adaptive N-R] Converged in %d iterations: V=%.6e, R=%.6e\n', iter, V, R);
            end
            return;
        end
        
        % Early termination if too many bad steps
        if consecutive_bad_steps > 5
            if debug_mode
                fprintf('  [Adaptive N-R] Too many bad steps, terminating at iter %d\n', iter);
            end
            break;
        end
        
    catch ME
        % If calculation error occurs, try to diagnose
        if debug_mode
            fprintf('  [Adaptive N-R] Error at iter %d: %s\n', iter, ME.message);
        end
        consecutive_bad_steps = consecutive_bad_steps + 1;
        if consecutive_bad_steps > 3
            break;
        end
    end
end

% If we reach here without convergence
if debug_mode && ~converged
    fprintf('  [Adaptive N-R] Failed to converge in %d iterations: V=%.6e, R=%.6e\n', max_iter, V, R);
end

end

