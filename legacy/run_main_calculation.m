function run_main_calculation(Cpv, DHv, a, Gibbs_Tom, D_Al, m_Al, k0_Al, sigma, ...
    deltaT_min, deltaT_max, sampling_interval, split_point, ...
    lower_density_mult, upper_density_mult, low_end_density, ...
    high_end_density, C0_values, V_a, V_b, V_c, R_d, R_e, R_f, ...
    include_thermal_undercooling, use_analytical_approximation)
%% RUN_MAIN_CALCULATION - Execute comprehensive LGK dendritic growth calculations
% ========================================================================
%
% PURPOSE:
% This is the main computational engine for the LGK model, orchestrating
% multi-concentration dendritic growth calculations with parallel processing,
% progress monitoring, and comprehensive result generation.
%
% CALCULATION WORKFLOW:
% 1. Parameter validation and parallel pool initialization
% 2. DeltaT sequence generation with adaptive density control
% 3. Multi-concentration calculation (parallel or serial)
% 4. Iterative parameter refinement for convergence optimization
% 5. Result compilation and comprehensive visualization
% 6. Surface fitting and statistical analysis
%
% PHYSICAL IMPLEMENTATION:
% Solves the coupled LGK equations for each concentration:
% - Undercooling balance: ΔT = ΔTₜ + ΔTc + ΔTᵣ
% - Stability criterion: R = σΓ/(mGc - G)
% Using Newton-Raphson iteration with adaptive initial guess refinement
%
% PARALLEL PROCESSING STRATEGY:
% - Automatic detection of optimal worker pool size
% - Load balancing across concentration values
% - Progress synchronization with GUI feedback
% - Graceful fallback to serial computation
%
% INPUTS:
%   Material Parameters:
%     Cpv - Volumetric specific heat [J/m³·K]
%     DHv - Volumetric latent heat [J/m³]
%     a - Thermal diffusivity [m²/s]
%     Gibbs_Tom - Gibbs-Thomson coefficient [K·m]
%     D_Al - Interdiffusion coefficient [m²/s]
%     m_Al - Liquidus slope [K/wt%]
%     k0_Al - Partition coefficient [-]
%     sigma - Stability constant [-]
%
%   Calculation Parameters:
%     deltaT_min/max - Undercooling range [K]
%     sampling_interval - Base spacing for deltaT sequence [K]
%     split_point - Transition point for two-region sampling [K]
%     density_mult - Region-specific sampling density multipliers [-]
%     C0_values - Array of solute concentrations [wt%]
%
%   Initial Guess Parameters:
%     V_a, V_b, V_c - Velocity function coefficients
%     R_d, R_e, R_f - Radius function coefficients
%     include_thermal_undercooling - Include thermal effects [logical]
%     use_analytical_approximation - Use theoretical vs parametric guess [logical]
%
% OUTPUTS:
%   Results stored in base workspace:
%     export_data - Comprehensive results structure for Excel export
%     all_C0_final_results - Final V,R values for each concentration
%     Material parameters and calculation metadata
%
% DEPENDENCIES:
%   calculate_VR.m - Core Newton-Raphson solver
%   create_deltaT_sequence.m - Adaptive undercooling sequence generation
%   create_analytical_guess_functions.m - Theoretical initial guess functions
%   finalize_plots.m - Comprehensive visualization and surface fitting
% ========================================================================
tic

% Add analytical approximation functions to parallel workers
if use_analytical_approximation
    % Ensure analytical functions are available in parallel workers  
    current_path = pwd;
    if ~isempty(gcp('nocreate'))
        pctRunOnAll(['addpath(''' current_path ''')']);
        fprintf('Added current path to parallel workers for analytical functions\n');
    end
    
    % Also ensure the function is available in the main thread
    if exist('create_analytical_guess_functions', 'file') ~= 2
        error('analytical_approximation_functions.m file not found in current directory: %s', current_path);
    end
end


% Ensure that output is not suppressed during normal calculation
global SUPPRESS_OUTPUT;
SUPPRESS_OUTPUT = false;

% Initialize progress tracking
global STOP_CALCULATION;
if isempty(STOP_CALCULATION)
    STOP_CALCULATION = false;
end

% Get GUI handles for progress updates
fig = findobj('Tag', 'MaterialControlPanel');
if isempty(fig)
    fig = findobj('Name', 'Material Parameter Control Panel');
end
gui_handles = [];
if ~isempty(fig)
    gui_handles = get(fig, 'UserData');
end

% Calculate total points FIRST
num_C0 = length(C0_values);

% Check if SCN material is selected (based on parameter values)
is_SCN_material = (deltaT_min == 0.5 && deltaT_max == 0.9 && sampling_interval == 0.4);

if is_SCN_material
    % For SCN material, use only two deltaT values
    deltaT_sequence = [0.5, 0.9];
    fprintf('SCN material detected: Using fixed deltaT values [0.5, 0.9]\n');
else
    % For other materials, use the standard deltaT sequence
    deltaT_sequence = create_deltaT_sequence('deltaT_min', deltaT_min, 'deltaT_max', deltaT_max, ...
        'sampling_interval', sampling_interval, 'split_point', split_point, ...
        'lower_density_mult', lower_density_mult, 'upper_density_mult', upper_density_mult, ...
        'low_end_density', low_end_density, 'high_end_density', high_end_density);
end

points_per_c0 = length(deltaT_sequence);
total_points = num_C0 * points_per_c0;

% Get material name from GUI
fig = findobj('Tag', 'MaterialControlPanel');
if isempty(fig)
    fig = findobj('Name', 'Material Parameter Control Panel');
end
material_name = 'Unknown';
if ~isempty(fig)
    gui_handles_temp = get(fig, 'UserData');
    if isstruct(gui_handles_temp) && isfield(gui_handles_temp, 'material_popup')
        material_names = {'AZ91', 'Al-4wt%Cu', 'Al-Cu', 'Al-Fe', 'Sn-Ag', 'Sn-Cu', 'SCN-Acetone(LGK 1984)', 'Mg-Alloy(Lin 2009)', 'Custom'};
        material_selection = get(gui_handles_temp.material_popup, 'Value');
        material_name = material_names{material_selection};
        % Replace special characters in the string to yield a valid field name
        material_name = strrep(material_name, '-', '_');
        material_name = strrep(material_name, '%', 'pct');
        material_name = strrep(material_name, '(', '_');
        material_name = strrep(material_name, ')', '_');
        material_name = strrep(material_name, ' ', '_');
    end
end

% Initialize progress manager with parallel pool status
% Pool should be ready since it was initialized in GUI callback
pool = gcp('nocreate');
pool_already_ready = ~isempty(pool);
if pool_already_ready
    fprintf('Progress manager: Parallel pool ready (%d workers)\n', pool.NumWorkers);
else
    fprintf('Progress manager: No parallel pool - will use serial computation\n');
end
progress_manager = init_progress_manager(material_name, total_points, gui_handles, pool_already_ready);

% Check if parallel pool already exists to adjust time estimates
pool = gcp('nocreate');
if ~isempty(pool)
    actual_workers = pool.NumWorkers;
    parallel_setup_time = 0; % No setup time needed - pool already ready
    pool_status = 'ready';
    fprintf('Parallel pool already available (%d workers) - no initialization time needed\n', actual_workers);
else
    assumed_workers = min(8, feature('numcores'));
    actual_workers = assumed_workers;
    parallel_setup_time = 0.5; % Reduced time since pool may be initializing in background
    pool_status = 'initializing or unavailable';
    fprintf('Parallel pool not ready - may need initialization time\n');
end

time_per_point = 1 / actual_workers; % seconds per point
estimated_total_time = parallel_setup_time + (total_points * time_per_point) / 60; % in minutes

fprintf('Pre-calculation estimates:\n');
fprintf('  Total points: %d\n', total_points);
fprintf('  Available parallel workers: %d\n', actual_workers);
fprintf('  Pool status: %s\n', pool_status);
fprintf('  Estimated setup time: %.1f minutes\n', parallel_setup_time);
fprintf('  Estimated total time: %.1f minutes\n', estimated_total_time);



% Record time before checking parallel status
parallel_check_start = now;
fprintf('\nChecking parallel computing status...\n');

% Check if parallel pool is ready (should be ready since initialized in GUI)
pool_was_ready = ~isempty(gcp('nocreate'));
if pool_was_ready
    fprintf('✓ Parallel pool is ready - proceeding with calculation\n');
else
    fprintf('⚠ No parallel pool found - will use serial computation\n');
end

% Record time after parallel check - use consistent variable name
parallel_check_end = now;
parallel_init_time = (parallel_check_end - parallel_check_start) * 24 * 60; % in minutes

if pool_was_ready
    fprintf('Parallel pool was already ready - no wait time (%.3f minutes)\n', parallel_init_time);
else
    fprintf('Parallel pool setup/wait completed in %.2f minutes\n', parallel_init_time);
end

% Update progress after parallel initialization
% Store the parallel wait time in progress manager
progress_manager.parallel_init_time = parallel_init_time;  % Keep the original field name for compatibility
progress_manager = update_progress_checkpoint(progress_manager, 'parallel_init_complete', 0, 'Parallel Computing Initialized');

% Get actual pool information
pool = gcp('nocreate');
if ~isempty(pool)
    actual_workers = pool.NumWorkers;
else
    actual_workers = 1; % Serial computation
end

% Update time estimate with actual workers
actual_time_per_point = 1 / (actual_workers); % seconds per point
actual_estimated_time = parallel_init_time + (total_points * actual_time_per_point) / 60; % in minutes

fprintf('Updated estimates with %d workers:\n', actual_workers);
fprintf('  Time per point: %.2f seconds\n', actual_time_per_point);
fprintf('  Parallel wait time: %.2f minutes\n', parallel_init_time);  
fprintf('  Total estimated time: %.1f minutes\n', actual_estimated_time);


% Store original material parameters
material_params = struct();
material_params.Cpv = Cpv;
material_params.DHv = DHv;
material_params.a = a;  % Thermal diffusivity
material_params.Gibbs_Tom = Gibbs_Tom;
material_params.D_Al = D_Al;
material_params.m_Al = m_Al;
material_params.k0_Al = k0_Al;
material_params.sigma = sigma;
material_params.include_thermal_undercooling = include_thermal_undercooling;



fprintf('\nStarting calculation with selected material parameters...\n');
fprintf('Material Parameters:\n');
fprintf('  Cpv = %.6e J/m³·K\n', Cpv);
fprintf('  DHv = %.6e J/m³\n', DHv);
fprintf('  a = %.6e m²/s\n', material_params.a);
fprintf('  Gibbs-Thomson = %.6e K·m\n', Gibbs_Tom);
fprintf('  D_Al = %.6e m²/s\n', D_Al);
fprintf('  m_Al = %.6f K/wt%%\n', m_Al);
fprintf('  k0_Al = %.6f\n', k0_Al);
fprintf('  σ = %.6f\n', sigma);
fprintf('V_guess Function Parameters:\n');
if use_analytical_approximation
    fprintf('Using analytical approximation (Eqs. 8.91 & 8.92) - Initial parametric values: V_a=%.6f, V_b=%.6f, V_c=%.6f\n', V_a, V_b, V_c);
    fprintf('R_guess parametric values: R_d=%.6e, R_e=%.6f, R_f=%.6f\n', R_d, R_e, R_f);
else
    fprintf('  V = (ΔT/(%.1f×C0^%.3f))^%.3f\n', V_a, V_b, V_c);
end

% Maximum iterations
max_iterations = 10;

% Store all results for different C0 values
all_C0_results = cell(num_C0, 1);
all_C0_final_results = cell(num_C0, 1);
all_C0_min_cnvg_DT = cell(num_C0, 1);

% Define colors for different C0 values
colors = lines(num_C0);

fprintf('\nStarting multi-concentration calculation...\n');
fprintf('C0 values: %.1f to %.1f\n', C0_values(1), C0_values(end));

fprintf('All C0 values: ');
for i = 1:length(C0_values)
    fprintf('%.3f ', C0_values(i));
    if mod(i, 10) == 0
        fprintf('\n               ');
    end
end
fprintf('\n');

if use_analytical_approximation
    fprintf('Using analytical approximation (Eqs. 8.91 & 8.92)\n');
else
    fprintf('Using parametric guess functions\n');
end

% Initialize global R_guess parameters that will be updated across C0 values
global_R_d = R_d;
global_R_e = R_e;
global_R_f = R_f;

% Initialize actual completed points counter
actual_completed_points = 0;

% Main loop continues from here...

% Main loop for different C0 values
% Determine if parallel processing is beneficial
pool = gcp('nocreate');
if ~isempty(pool) && num_C0 >= 2
    fprintf('Using parallel processing for %d C0 values with %d workers\n', num_C0, pool.NumWorkers);
    % ... (rest of the parallel processing code)

    % Parallel processing of C0 values
    all_C0_results_par = cell(num_C0, 1);
    all_C0_final_results_par = cell(num_C0, 1);
    all_C0_min_cnvg_DT_par = cell(num_C0, 1);

    % Convert variables for parallel workers
    par_material_params = material_params;
    par_deltaT_sequence = deltaT_sequence;
    par_max_iterations = max_iterations;
    par_global_R_d = global_R_d;
    par_global_R_e = global_R_e;
    par_global_R_f = global_R_f;
    par_use_analytical = use_analytical_approximation;

    fprintf('Starting parallel computation across C0 values...\n');
    
    parfor c0_idx = 1:num_C0
    C0 = C0_values(c0_idx);
    fprintf('Worker processing C0 = %.1f\n', C0);
    
    [worker_results, worker_final, worker_min_dt] = process_single_c0_worker(...
            c0_idx, C0, par_material_params, par_deltaT_sequence, par_max_iterations, ...
            V_a, V_b, V_c, par_global_R_d, par_global_R_e, par_global_R_f, par_use_analytical);
    
        all_C0_results_par{c0_idx} = worker_results;
        all_C0_final_results_par{c0_idx} = worker_final;
        all_C0_min_cnvg_DT_par{c0_idx} = worker_min_dt;
    
        
        fprintf('C0 = %.1f completed by worker\n', C0);
    end
        

    
    fprintf('Parallel computation completed for all C0 values\n');

    % Copy parallel results back
    all_C0_results = all_C0_results_par;
    all_C0_final_results = all_C0_final_results_par;
    all_C0_min_cnvg_DT = all_C0_min_cnvg_DT_par;

    fprintf('Parallel computation completed for all C0 values\n');

    % Update progress for parallel completion
    progress_manager = update_progress_checkpoint(progress_manager, 'calculation_progress', ...
        total_points, 'All C0 Values Complete');

else
    % Original serial processing
    fprintf('Using serial processing for C0 values\n');

    % Main loop for different C0 values
    for c0_idx = 1:num_C0
        % Check for stop signal
        if STOP_CALCULATION
            fprintf('Calculation stopped by user at C0 = %.1f\n', C0_values(c0_idx));
            break;
        end
        C0 = C0_values(c0_idx);
        current_color = colors(c0_idx, :);

        fprintf('\n========================================\n');
        fprintf('Processing C0 = %.1f (%.0f/%.0f)\n', C0, c0_idx, num_C0);
        fprintf('========================================\n');

        % Reset initial guess parameters for each C0
        V_a_current = V_a;
        V_b_current = V_b;
        V_c_current = V_c;
        % Use global R parameters (updated from previous C0)
        R_d_current = global_R_d;
        R_e_current = global_R_e;
        R_f_current = global_R_f;

        fprintf('Using R_guess parameters from previous C0: R_d = %.6e, R_e = %.6f, R_f = %.6f\n', global_R_d, global_R_e, global_R_f);

        % Store iteration results for current C0
        all_iterations = cell(max_iterations, 1);
        convergence_history = zeros(max_iterations, 1);
        min_cnvg_DT = zeros(max_iterations, 1);

        % Try different V_a scaling factors if initial attempt fails
        scaling_factors = [1, 2, 0.5, 4, 0.25, 8, 0.125, 16, 0.0625, 100, 1e-2, 1000, 1e-3];
        best_result = [];
        best_scaling = 1;

        % Store the best test result and its parameters
        best_test_result = [];
        best_converged_count = 0;

        for scale_idx = 1:length(scaling_factors)
            scale_factor = scaling_factors(scale_idx);

            if scale_idx > 1
                fprintf('\n--- Trying V_a scaling factor: %.4f ---\n', scale_factor);
            end

            % Reset parameters with scaling
            V_a_test = V_a * scale_factor;
            V_b_test = V_b;
            V_c_test = V_c;
            R_d_test = global_R_d;
            R_e_test = global_R_e;
            R_f_test = global_R_f;

            % Quick test with scaling factor - only test first few points
            V_guess_func_test = @(deltaT) (deltaT./(V_a_test * C0.^V_b_test)).^V_c_test;
            R_guess_func_test = @(deltaT) R_d_test * (C0.^R_e_test) .* (deltaT.^R_f_test);

            % Use only first 5 deltaT values for quick testing
            test_deltaT_sequence = deltaT_sequence(1:min(5, length(deltaT_sequence)));
            fprintf('Testing with deltaT values: ');
            for j = 1:length(test_deltaT_sequence)
                fprintf('%.2f ', test_deltaT_sequence(j));
            end
            fprintf('\n');

            [test_result, ~, ~] = calculate_VR(V_guess_func_test, R_guess_func_test, C0, test_deltaT_sequence, material_params);

            % Update progress after testing phase (Calibration Point 1)
            if ~isempty(test_result)
                test_points = size(test_result, 1);
                actual_completed_points = actual_completed_points + test_points;
                
                
            end

            % Store this test result if it's better than previous ones
            if size(test_result, 1) > best_converged_count
                best_test_result = test_result;
                best_converged_count = size(test_result, 1);
                best_V_a = V_a_test;
                best_V_b = V_b_test;
                best_V_c = V_c_test;
                best_R_d = R_d_test;
                best_R_e = R_e_test;
                best_R_f = R_f_test;
                best_scaling = scale_factor;
            end

            % Check if current result is good enough to proceed
            if size(test_result, 1) >= 3
                fprintf('Scaling factor %.4f successful, proceeding with formal calculation\n', scale_factor);
                break;
            end

            % Handle case where no scaling factor works
            if scale_idx == length(scaling_factors)
                if best_converged_count == 0
                    fprintf('No scaling factor worked, skipping C0 = %.1f\n', C0);
                    all_C0_results{c0_idx} = {[]};
                    all_C0_final_results{c0_idx} = [];
                    all_C0_min_cnvg_DT{c0_idx} = [];
                    continue;  % Skip to next C0
                else
                    fprintf('Using best available result with scaling %.4f (%d points)\n', best_scaling, best_converged_count);
                    % Use the best result found, even if insufficient for iteration
                    all_C0_results{c0_idx} = {best_test_result};
                    all_C0_final_results{c0_idx} = best_test_result;
                    all_C0_min_cnvg_DT{c0_idx} = best_test_result(end, 1);
                    fprintf('\nC0 = %.1f completed with limited results. Final converged points: %d\n', C0, best_converged_count);
                    continue;  % Skip to next C0
                end
            end
        end

        % Check if we have a successful test result to inherit
        if best_converged_count >= 3
            % Use the best test result parameters
            V_a_current = best_V_a;
            V_b_current = best_V_b;
            V_c_current = best_V_c;
            R_d_current = best_R_d;
            R_e_current = best_R_e;
            R_f_current = best_R_f;

            fprintf('\n--- Formal Calculation Starting ---\n');
            fprintf('Inheriting %d converged test points (deltaT: ', best_converged_count);
            for j = 1:best_converged_count
                fprintf('%.2f', best_test_result(j, 1));
                if j < best_converged_count
                    fprintf(', ');
                end
            end
            fprintf(')\n');

            if best_converged_count < length(deltaT_sequence)
                fprintf('Continuing calculation from deltaT = %.2f onwards\n', deltaT_sequence(best_converged_count + 1));
                % Determine which deltaT points still need to be calculated
                remaining_deltaT_sequence = deltaT_sequence(best_converged_count + 1:end);
            else
                fprintf('All deltaT points were calculated during testing\n');
                remaining_deltaT_sequence = [];
            end

        else
            fprintf('\nNo suitable test results found, skipping C0 = %.1f\n', C0);
            all_C0_results{c0_idx} = {[]};
            all_C0_final_results{c0_idx} = [];
            all_C0_min_cnvg_DT{c0_idx} = [];
            continue;  % Skip to next C0
        end


        % Iterative calculation for current C0
        for iter = 1:max_iterations
            % Check for stop signal
            if STOP_CALCULATION
                fprintf('Calculation stopped by user during iteration %d for C0 = %.1f\n', iter, C0);
                break;
            end
            fprintf('\n=== C0=%.1f, Iteration %d (V_a scaling: %.4f) ===\n', C0, iter, best_scaling);

            % Construct current guess functions
            % Check if using analytical approximation
            if use_analytical_approximation
                % Use analytical approximation based on Eqs. (8.91) and (8.92)
                [V_guess_func, R_guess_func] = create_analytical_guess_functions(C0, Gibbs_Tom, m_Al, k0_Al, D_Al);
                fprintf('Using analytical approximation (Eqs. 8.91 & 8.92) for C0=%.1f\n', C0);
            elseif is_SCN_material
                % Use fitted polynomial and power function interpolation for SCN-Acetone
                V_guess_func = @(deltaT) arrayfun(@(dt) get_SCN_fitted_V_guess(C0, dt), deltaT);
                R_guess_func = @(deltaT) arrayfun(@(dt) get_SCN_fitted_R_guess(C0, dt), deltaT);
                
                fprintf('Using polynomial interpolation for SCN-Acetone material (C0=%.3f)\n', C0);
            else
                % Use original parametric form for other materials
                V_guess_func = @(deltaT) (deltaT./(V_a_current * C0.^V_b_current)).^V_c_current;
                R_guess_func = @(deltaT) R_d_current * (C0.^R_e_current) .* (deltaT.^R_f_current);
            end
            fprintf('Current guess function parameters:\n');
            fprintf('V_guess: (ΔT/(%.6f * C0^%.6f))^%.6f\n', V_a_current, V_b_current, V_c_current);
            fprintf('R_guess: %.6e * C0^%.6f * ΔT^%.6f\n', R_d_current, R_e_current, R_f_current);
            fprintf('For C0 = %.1f: V_guess = (ΔT/%.6f)^%.6f\n', C0, V_a_current * C0^V_b_current, V_c_current);
            fprintf('Calculating deltaT values: ');

            %         % Call calculation function
            %         [result, converged, last_diverged_deltaT] = calculate_VR(V_guess_func, R_guess_func, C0, deltaT_sequence, material_params);
            % Calculate only the remaining deltaT points
            if ~isempty(remaining_deltaT_sequence)
                [remaining_result, converged, last_diverged_deltaT] = calculate_VR(V_guess_func, R_guess_func, C0, remaining_deltaT_sequence, material_params);

                % Combine inherited test results with new calculations
                if ~isempty(remaining_result)
                    result = [best_test_result; remaining_result];
                    fprintf('\nCombined results: %d inherited + %d newly calculated = %d total points\n', ...
                        size(best_test_result, 1), size(remaining_result, 1), size(result, 1));
                else
                    result = best_test_result;
                    fprintf('\nUsing only inherited test results (%d points)\n', size(result, 1));
                end
            else
                % All points were already calculated in testing
                result = best_test_result;
                fprintf('\nAll required points were calculated during testing (%d points)\n', size(result, 1));
                converged = true;
                last_diverged_deltaT = NaN;
            end

            % Calculate only the remaining deltaT points
            if ~isempty(remaining_deltaT_sequence)
                [remaining_result, converged, last_diverged_deltaT] = calculate_VR(V_guess_func, R_guess_func, C0, remaining_deltaT_sequence, material_params);

                % Combine inherited test results with new calculations
                if ~isempty(remaining_result)
                    result = [best_test_result; remaining_result];
                    fprintf('\nCombined results: %d inherited + %d newly calculated = %d total points\n', ...
                        size(best_test_result, 1), size(remaining_result, 1), size(result, 1));
                else
                    result = best_test_result;
                    fprintf('\nUsing only inherited test results (%d points)\n', size(result, 1));
                end
            else
                % All points were already calculated in testing
                result = best_test_result;
                fprintf('\nAll required points were calculated during testing (%d points)\n', size(result, 1));
                converged = true;
                last_diverged_deltaT = NaN;
            end

            % Update progress after iteration completion (Calibration Point 2)
            if ~isempty(result)
                current_iteration_points = size(result, 1);
                % Only count new points, subtract previously counted test points if they were included
                if iter == 1
                    new_points = current_iteration_points - size(best_test_result, 1);
                else
                    new_points = current_iteration_points - size(all_iterations{iter-1}, 1);
                end
                actual_completed_points = actual_completed_points + max(0, new_points);
                

                
            end

            % Store results
            all_iterations{iter} = result;
            convergence_history(iter) = size(result, 1);

            % Record minimum convergent deltaT
            if ~isempty(result)
                min_cnvg_DT(iter) = result(end, 1);  % Last converged deltaT
                fprintf('\nIteration %d completed: Minimum convergent deltaT = %.2f\n', iter, min_cnvg_DT(iter));
            else
                min_cnvg_DT(iter) = NaN;
                fprintf('\nIteration %d: No convergent points\n', iter);
            end

            fprintf('Converged points in this iteration: %d\n', size(result, 1));
            if ~isnan(last_diverged_deltaT)
                fprintf('Last diverged deltaT: %.2f\n', last_diverged_deltaT);
            else
                fprintf('All deltaT converged\n');
            end

            % Check termination condition 1: All converged
            if size(result, 1) == length(deltaT_sequence)
                fprintf('\n*** Termination condition 1 reached: All deltaT converged ***\n');
                break;
            end

            % Check termination condition 2: No improvement in past two iterations
            if iter >= 3
                if convergence_history(iter) <= convergence_history(iter-1) && ...
                        convergence_history(iter-1) <= convergence_history(iter-2)
                    fprintf('\n*** Termination condition 2 reached: No improvement in convergence ***\n');
                    break;
                end
            end

            % If too few points for fitting, exit
            if size(result, 1) < 3
                fprintf('\n*** Too few converged points for fitting, exiting iteration ***\n');
                break;
            end

            % Update inherited results for next iteration - keep the better results
            if size(result, 1) > size(best_test_result, 1)
                best_test_result = result;
                best_converged_count = size(result, 1);
                if best_converged_count < length(deltaT_sequence)
                    remaining_deltaT_sequence = deltaT_sequence(best_converged_count + 1:end);
                else
                    remaining_deltaT_sequence = [];
                end
            end

            % Fit new guess function parameters
            fprintf('Fitting new guess functions...\n');

            % R_guess_func fitting: R = d * C0^e * deltaT^f
            % Taking log: log(R) = log(d) + e*log(C0) + f*log(deltaT)
            deltaT_data = result(:, 1);
            R_data = result(:, 3);

            log_deltaT = log(deltaT_data);
            log_C0 = log(repmat(C0, length(deltaT_data), 1));
            log_R = log(R_data);

            % Create design matrix for multiple linear regression
            X_R = [ones(length(log_R), 1), log_C0, log_deltaT];

            % Solve linear system: log(R) = X_R * coeffs
            R_coeffs = X_R \ log_R;

            % Extract parameters
            R_d_new = exp(R_coeffs(1));  % d parameter
            R_e_new = R_coeffs(2);       % e parameter
            R_f_new = R_coeffs(3);       % f parameter

            % V_guess_func fitting: use last 50% data points
            % NEW form: V = (deltaT/(a*C0^b))^c
            % Rearranging: deltaT = a*C0^b * V^(1/c)
            % For fixed C0: deltaT = K * V^(1/c), where K = a*C0^b
            % So: V = (deltaT/K)^c
            % Taking log: log(V) = c*log(deltaT) - c*log(K)
            % We fit: log(V) = slope*log(deltaT) + intercept
            % where slope = c and intercept = -c*log(K) = -c*log(a) - c*b*log(C0)

            n_points = size(result, 1);
            start_idx = max(1, ceil(n_points * 0.5));

            deltaT_data_V = result(start_idx:end, 1);
            V_data = result(start_idx:end, 2);

            log_deltaT_V = log(deltaT_data_V);
            log_V = log(V_data);

            % Fit log(V) vs log(deltaT)
            V_coeffs = polyfit(log_deltaT_V, log_V, 1);
            V_c_new = V_coeffs(1);  % This is the new c parameter
            intercept = V_coeffs(2);

            % From intercept = -c*log(a) - c*b*log(C0), solve for a:
            % log(a) = -(intercept + c*b*log(C0))/c
            V_a_new = exp(-(intercept + V_c_new * V_b_current * log(C0)) / V_c_new);

            % Keep V_b unchanged for stability, or allow slight adjustment
            V_b_new = V_b_current;

            fprintf('New fitted parameters:\n');
            fprintf('V_guess: (ΔT/(%.6f * C0^%.6f))^%.6f\n', V_a_new, V_b_new, V_c_new);
            fprintf('R_guess: %.6e * C0^%.6f * ΔT^%.6f\n', R_d_new, R_e_new, R_f_new);

            % Update parameters
            V_a_current = V_a_new;
            V_b_current = V_b_new;
            V_c_current = V_c_new;
            R_d_current = R_d_new;
            R_e_current = R_e_new;
            R_f_current = R_f_new;
        end

        % Store results for current C0
        all_C0_results{c0_idx} = all_iterations;
        all_C0_final_results{c0_idx} = all_iterations{iter};
        all_C0_min_cnvg_DT{c0_idx} = min_cnvg_DT(1:iter);

        % Plot results for current C0
        final_result = all_iterations{iter};

        %     if ~isempty(final_result)
        %         plot_current_results(final_result, current_color, C0, iter);
        %     end

        fprintf('\nC0 = %.1f completed. Final converged points: %d\n', C0, size(final_result, 1));
        
        % Update progress after each C0 completion in serial mode
        completed_points_so_far = c0_idx * points_per_c0;
        progress_manager = update_progress_checkpoint(progress_manager, 'calculation_progress', ...
            completed_points_so_far, sprintf('C0=%.1f Complete (%d/%d)', C0, c0_idx, num_C0));
        

        % Store results for current C0
        all_C0_results{c0_idx} = all_iterations;
        all_C0_final_results{c0_idx} = all_iterations{iter};
        all_C0_min_cnvg_DT{c0_idx} = min_cnvg_DT(1:iter);

        % Update global R parameters for next C0 (NEW CODE)
        if ~isempty(all_iterations{iter}) && size(all_iterations{iter}, 1) >= 3
            % Use the final optimized R parameters from this C0
            global_R_d = R_d_current;
            global_R_e = R_e_current;
            global_R_f = R_f_current;
            fprintf('Updated global R parameters: R_d = %.6e, R_e = %.6f, R_f = %.6f\n', global_R_d, global_R_e, global_R_f);
        else
            fprintf('Keeping previous R parameters due to insufficient data\n');
        end
    end % End of serial C0 loop
end % End of parallel/serial choice



% Display timing information
timer = toc;
fprintf('\nTotal computation time: %.2f seconds\n', timer);

% Store all results in a comprehensive structure for plotting and Excel export
export_data = struct();
export_data.C0_values = C0_values;
export_data.all_results = {};
export_data.material_name = material_name;  % Use the correct material_name variable
export_data.material_params = material_params;
export_data.calculation_time = timer;
export_data.deltaT_params = struct('deltaT_min', deltaT_min, 'deltaT_max', deltaT_max, ...
    'sampling_interval', sampling_interval, 'split_point', split_point, ...
    'lower_density_mult', lower_density_mult, 'upper_density_mult', upper_density_mult, ...
    'low_end_density', low_end_density, 'high_end_density', high_end_density);

% Store INITIAL guess parameters (CRITICAL for showing guess functions in plots)
export_data.initial_guess_params = struct('V_a', V_a, 'V_b', V_b, 'V_c', V_c, ...
    'R_d', R_d, 'R_e', R_e, 'R_f', R_f);

% Organize calculation results for each C0
for c0_idx = 1:length(C0_values)
    if ~isempty(all_C0_final_results{c0_idx})
        export_data.all_results{c0_idx} = struct();
        export_data.all_results{c0_idx}.C0 = C0_values(c0_idx);
        export_data.all_results{c0_idx}.deltaT = all_C0_final_results{c0_idx}(:, 1);
        export_data.all_results{c0_idx}.R = all_C0_final_results{c0_idx}(:, 3);
        export_data.all_results{c0_idx}.V = all_C0_final_results{c0_idx}(:, 2);
    else
        export_data.all_results{c0_idx} = struct();
        export_data.all_results{c0_idx}.C0 = C0_values(c0_idx);
        export_data.all_results{c0_idx}.deltaT = [];
        export_data.all_results{c0_idx}.R = [];
        export_data.all_results{c0_idx}.V = [];
    end
end

% Store all data in base workspace for GUI access and plotting
assignin('base', 'export_data', export_data);
assignin('base', 'all_C0_final_results', all_C0_final_results);
assignin('base', 'C0_values', C0_values);
assignin('base', 'material_params', material_params);
assignin('base', 'calculation_time', timer);
assignin('base', 'current_material_name', material_name);

fprintf('Export_data with initial parameters saved: V_a=%.3g, V_b=%.3g, V_c=%.3g, R_d=%.3g, R_e=%.3g, R_f=%.3g\n', ...
        V_a, V_b, V_c, R_d, R_e, R_f);
fprintf('Calculation results stored for Excel export.\n');

% Finalize all plots (now export_data is properly available)
finalize_plots(all_C0_final_results, C0_values, colors, material_name);

% Adjust figure layout
adjust_figure_layout();



% Final progress update using new progress manager
progress_manager = update_progress_checkpoint(progress_manager, 'calculation_complete', total_points, 'Calculation Complete');

% Calculate and display final timing
final_time = now;
total_elapsed_time = (final_time - progress_manager.start_time) * 24 * 60;
fprintf('Total elapsed time: %.2f minutes\n', total_elapsed_time);
if progress_manager.pool_already_ready
    fprintf('  ✓ Parallel pool was pre-initialized - saved setup time\n');
    fprintf('  📊 Pure calculation time: %.2f minutes\n', total_elapsed_time);
else
    if isfield(progress_manager, 'parallel_init_time') && progress_manager.parallel_init_time > 0
        fprintf('  ⏱ Parallel pool setup time: %.2f minutes\n', progress_manager.parallel_init_time);
        fprintf('  📊 Pure calculation time: %.2f minutes\n', total_elapsed_time - progress_manager.parallel_init_time);
    else
        fprintf('  📊 Calculation time: %.2f minutes\n', total_elapsed_time);
    end
end

% Keep parallel pool running for next calculation
pool = gcp('nocreate');
if ~isempty(pool)
    fprintf('Parallel pool kept running (%d workers) for next calculation\n', pool.NumWorkers);
else
    fprintf('No parallel pool to maintain\n');
end

% Clean up temporary progress files
temp_files = dir('temp_progress_*.mat');
for i = 1:length(temp_files)
    try
        delete(temp_files(i).name);
    catch
        % Silent cleanup
    end
end
if ~isempty(temp_files)
    fprintf('Cleaned up %d temporary progress files\n', length(temp_files));
end

% At the end of the function, restore default settings
SUPPRESS_OUTPUT = false;

end

function finalize_plots(all_C0_final_results, C0_values, colors, material_name)
%% FINALIZE_PLOTS - Generate comprehensive visualization suite
%
% PURPOSE:
% Creates complete visualization package including 2D plots with initial
% guess function overlays, concentration-based color mapping, and 
% comprehensive figure layout management.
%
% VISUALIZATION STRATEGY:
% - Concentration-based jet colormap for consistent data representation
% - Initial guess function overlays for model validation
% - Multiple scaling options (linear/logarithmic) for different data ranges
% - Automated figure positioning for optimal screen utilization
%
% PLOT TYPES GENERATED:
% Figure 1: V vs C₀ (Linear) with initial guess overlays
% Figure 2: R vs C₀ (Linear) with initial guess overlays  
% Figure 3: V vs C₀ (Log) for wide dynamic range visualization
% Figures 6-7: ΔT-based plots with concentration differentiation
%
% COLOR MAPPING ALGORITHM:
% Uses jet colormap with concentration-proportional indexing:
% - C₀_normalized = (C₀ - C₀_min)/(C₀_max - C₀_min)
% - Color_index = round(C₀_normalized × 255) + 1
% Ensures consistent color representation across all plots

% Generate a colour map based on concentration values – map the actual concentration range to the jet colour scheme
num_C0 = length(C0_values);
if num_C0 > 1
    % Obtain the minimum and maximum concentration values
    C0_min = min(C0_values);
    C0_max = max(C0_values);

    % Calculate the relative position (between 0 and 1) of each concentration within the range
    C0_normalized = (C0_values - C0_min) / (C0_max - C0_min);

    % Retrieve the full jet colour map
    jet_colormap = jet(256);

    % Interpolate the corresponding colour based on the normalised concentration position
    concentration_colors = zeros(num_C0, 3);
    for i = 1:num_C0
        % Map the normalised position to a colour map index (from 1 to 256)
        color_index = round(C0_normalized(i) * 255) + 1;
        color_index = max(1, min(256, color_index)); % Ensure the index is within the valid range
        concentration_colors(i, :) = jet_colormap(color_index, :);
    end
else
    % Use blue for a single concentration value
    concentration_colors = [0, 0, 1];
end

% Obtain initial guess parameters for the function (multiple methods used to ensure success)
has_initial_params = false;
initial_V_a = NaN; initial_V_b = NaN; initial_V_c = NaN;
initial_R_d = NaN; initial_R_e = NaN; initial_R_f = NaN;

% Method 1: Try to get from base workspace export_data
try
    if evalin('base', 'exist(''export_data'', ''var'')')
        export_data = evalin('base', 'export_data');
        if isfield(export_data, 'initial_guess_params') && isstruct(export_data.initial_guess_params)
            igp = export_data.initial_guess_params;
            if isfield(igp, 'V_a') && isfield(igp, 'V_b') && isfield(igp, 'V_c') && ...
               isfield(igp, 'R_d') && isfield(igp, 'R_e') && isfield(igp, 'R_f')
                initial_V_a = igp.V_a; initial_V_b = igp.V_b; initial_V_c = igp.V_c;
                initial_R_d = igp.R_d; initial_R_e = igp.R_e; initial_R_f = igp.R_f;
                has_initial_params = true;
                fprintf('Method 1 SUCCESS: Got parameters from export_data\n');
            end
        end
    end
catch ME
    fprintf('Method 1 failed: %s\n', ME.message);
end

% Method 2: Try to get from GUI if Method 1 failed
if ~has_initial_params
    try
        fig = findobj('Tag', 'MaterialControlPanel');
        if ~isempty(fig)
            gui_handles = get(fig, 'UserData');
            if isstruct(gui_handles)
                initial_V_a = str2double(get(gui_handles.v_a_edit, 'String'));
                initial_V_b = str2double(get(gui_handles.v_b_edit, 'String'));
                initial_V_c = str2double(get(gui_handles.v_c_edit, 'String'));
                initial_R_d = str2double(get(gui_handles.rd_edit, 'String'));
                initial_R_e = str2double(get(gui_handles.re_edit, 'String'));
                initial_R_f = str2double(get(gui_handles.rf_edit, 'String'));
                
                % Check if we got valid numbers
                if ~isnan(initial_V_a) && ~isnan(initial_V_b) && ~isnan(initial_V_c) && ...
                   ~isnan(initial_R_d) && ~isnan(initial_R_e) && ~isnan(initial_R_f)
                    has_initial_params = true;
                    fprintf('Method 2 SUCCESS: Got parameters from GUI\n');
                end
            end
        end
    catch ME
        fprintf('Method 2 failed: %s\n', ME.message);
    end
end

% Method 3: Use default values if both methods failed
if ~has_initial_params
    initial_V_a = 52.8964; initial_V_b = 0.6; initial_V_c = 2.5;
    initial_R_d = 20.34e-6; initial_R_e = 0.25; initial_R_f = -1.25;
    has_initial_params = true;
    fprintf('Method 3 FALLBACK: Using default Al-Cu parameters\n');
end

fprintf('Final parameters: V_a=%.3g, V_b=%.3g, V_c=%.3g, R_d=%.3g, R_e=%.3g, R_f=%.3g\n', ...
        initial_V_a, initial_V_b, initial_V_c, initial_R_d, initial_R_e, initial_R_f);

% Retrieve all unique deltaT values for the legend
all_deltaT_values = [];
for c0_idx = 1:num_C0
    final_result = all_C0_final_results{c0_idx};
    if ~isempty(final_result) && size(final_result, 1) > 1
        all_deltaT_values = [all_deltaT_values; final_result(:,1)];
    end
end
unique_deltaT_values = unique(all_deltaT_values);
unique_deltaT_values = sort(unique_deltaT_values, 'descend'); % Sort in descending order

% Assign a colour to each deltaT value
deltaT_colors = jet(length(unique_deltaT_values));

% Figure 1: V vs C0 (Linear Scale) with guess functions
figure(1);
clf; hold on;

% Plot initial V guess function if available
if has_initial_params && ~isnan(initial_V_a) && ~isnan(initial_V_b) && ~isnan(initial_V_c)
    for deltaT_idx = 1:length(unique_deltaT_values)
        target_deltaT = unique_deltaT_values(deltaT_idx);
        current_color = deltaT_colors(deltaT_idx, :) * 0.5; % Darker for guess lines
        
        % Calculate initial guess function values
        C0_guess_range = linspace(min(C0_values), max(C0_values), 100);
        try
            V_initial_guess = (target_deltaT./(initial_V_a * C0_guess_range.^initial_V_b)).^initial_V_c;
            
            % Plot initial guess as dashed line
            plot(C0_guess_range, V_initial_guess, '--', 'Color', current_color, ...
                'LineWidth', 1.5, 'HandleVisibility', 'off');
        catch
            % Skip if calculation fails
        end
    end
    
    % Add legend entry for initial guess
    plot(NaN, NaN, '--', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 2, ...
        'DisplayName', sprintf('Initial: V=(ΔT/(%.2g×C₀^{%.3g}))^{%.3g}', initial_V_a, initial_V_b, initial_V_c));
end

% Plot actual data points
for deltaT_idx = 1:length(unique_deltaT_values)
    target_deltaT = unique_deltaT_values(deltaT_idx);
    current_color = deltaT_colors(deltaT_idx, :);
    
    % Collect all V values corresponding to each C0 at the given deltaT
    C0_data = [];
    V_data = [];
    
    for c0_idx = 1:num_C0
        final_result = all_C0_final_results{c0_idx};
        if ~isempty(final_result) && size(final_result, 1) > 1
            C0 = C0_values(c0_idx);
            
            % Find the data point closest to the target deltaT
            [~, closest_idx] = min(abs(final_result(:,1) - target_deltaT));
            if abs(final_result(closest_idx,1) - target_deltaT) < 0.1 % Allow a tolerance of 0.1 K
                C0_data = [C0_data; C0];
                V_data = [V_data; final_result(closest_idx, 2)];
            end
        end
    end
    
    % Plot the C0 vs V curve for this deltaT
    if length(C0_data) > 1
        [sorted_C0, sort_idx] = sort(C0_data);
        sorted_V = V_data(sort_idx);
        
        plot(sorted_C0, sorted_V, 'o', 'Color', current_color, ...
            'LineWidth', 1, 'MarkerSize', 5, 'DisplayName', sprintf('ΔT=%.1f K', target_deltaT));
        plot(sorted_C0, sorted_V, '-', 'Color', current_color, ...
            'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
end

xlabel('C₀ (wt%)');
ylabel('V (m s⁻¹)');
title(sprintf('%s - Velocity vs Concentration (Linear Scale)', material_name));
legend('Location', 'best');
grid on;

% Figure 2: R vs C0 (Linear Scale) with guess functions
figure(2);
clf; hold on;

% Plot initial R guess function if available
if has_initial_params && ~isnan(initial_R_d) && ~isnan(initial_R_e) && ~isnan(initial_R_f)
    for deltaT_idx = 1:length(unique_deltaT_values)
        target_deltaT = unique_deltaT_values(deltaT_idx);
        current_color = deltaT_colors(deltaT_idx, :) * 0.5; % Darker for guess lines
        
        % Calculate initial guess function values
        C0_guess_range = linspace(min(C0_values), max(C0_values), 100);
        try
            R_initial_guess = initial_R_d * (C0_guess_range.^initial_R_e) .* (target_deltaT.^initial_R_f);
            
            % Plot initial guess as dashed line
            plot(C0_guess_range, R_initial_guess, '--', 'Color', current_color, ...
                'LineWidth', 1.5, 'HandleVisibility', 'off');
        catch
            % Skip if calculation fails
        end
    end
    
    % Add legend entry for initial guess
    plot(NaN, NaN, '--', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 2, ...
        'DisplayName', sprintf('Initial: R=%.2g×C₀^{%.3g}×ΔT^{%.3g}', initial_R_d, initial_R_e, initial_R_f));
end

% Plot actual data points
for deltaT_idx = 1:length(unique_deltaT_values)
    target_deltaT = unique_deltaT_values(deltaT_idx);
    current_color = deltaT_colors(deltaT_idx, :);
    
    % Collect all R values corresponding to each C0 at the given deltaT
    C0_data = [];
    R_data = [];
    
    for c0_idx = 1:num_C0
        final_result = all_C0_final_results{c0_idx};
        if ~isempty(final_result) && size(final_result, 1) > 1
            C0 = C0_values(c0_idx);
            
            % Find the data point closest to the target deltaT
            [~, closest_idx] = min(abs(final_result(:,1) - target_deltaT));
            if abs(final_result(closest_idx,1) - target_deltaT) < 0.1 % Allow a tolerance of 0.1 K
                C0_data = [C0_data; C0];
                R_data = [R_data; final_result(closest_idx, 3)];
            end
        end
    end
    
    % Plot the C0 vs R curve for this deltaT
    if length(C0_data) > 1
        [sorted_C0, sort_idx] = sort(C0_data);
        sorted_R = R_data(sort_idx);
        
        plot(sorted_C0, sorted_R, 'o', 'Color', current_color, ...
            'LineWidth', 1, 'MarkerSize', 5, 'DisplayName', sprintf('ΔT=%.1f K', target_deltaT));
        plot(sorted_C0, sorted_R, '-', 'Color', current_color, ...
            'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
end

xlabel('C₀ (wt%)');
ylabel('R (m)');
title(sprintf('%s - Radius vs Concentration (Linear Scale)', material_name));
legend('Location', 'best');
grid on;

% Figure 3: V vs C0 (Log Scale) with guess functions
figure(3);
clf; hold on;

% Plot initial V guess function if available
if has_initial_params && ~isnan(initial_V_a) && ~isnan(initial_V_b) && ~isnan(initial_V_c)
    for deltaT_idx = 1:length(unique_deltaT_values)
        target_deltaT = unique_deltaT_values(deltaT_idx);
        current_color = deltaT_colors(deltaT_idx, :) * 0.5; % Darker for guess lines
        
        % Calculate initial guess function values
        C0_guess_range = linspace(min(C0_values), max(C0_values), 100);
        try
            V_initial_guess = (target_deltaT./(initial_V_a * C0_guess_range.^initial_V_b)).^initial_V_c;
            
            % Plot initial guess as dashed line
            plot(C0_guess_range, V_initial_guess, '--', 'Color', current_color, ...
                'LineWidth', 1.5, 'HandleVisibility', 'off');
        catch
            % Skip if calculation fails
        end
    end
    
    % Add legend entry for initial guess
    plot(NaN, NaN, '--', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 2, ...
        'DisplayName', sprintf('Initial: V=(ΔT/(%.2g×C₀^{%.3g}))^{%.3g}', initial_V_a, initial_V_b, initial_V_c));
end

% Plot actual data points
for deltaT_idx = 1:length(unique_deltaT_values)
    target_deltaT = unique_deltaT_values(deltaT_idx);
    current_color = deltaT_colors(deltaT_idx, :);
    
    % Collect all V values corresponding to each C0 at the given deltaT
    C0_data = [];
    V_data = [];
    
    for c0_idx = 1:num_C0
        final_result = all_C0_final_results{c0_idx};
        if ~isempty(final_result) && size(final_result, 1) > 1
            C0 = C0_values(c0_idx);
            
            % Find the data point closest to the target deltaT
            [~, closest_idx] = min(abs(final_result(:,1) - target_deltaT));
            if abs(final_result(closest_idx,1) - target_deltaT) < 0.1 % Allow a tolerance of 0.1 K
                C0_data = [C0_data; C0];
                V_data = [V_data; final_result(closest_idx, 2)];
            end
        end
    end
    
    % Plot the C0 vs V curve for this deltaT
    if length(C0_data) > 1
        [sorted_C0, sort_idx] = sort(C0_data);
        sorted_V = V_data(sort_idx);
        
        plot(sorted_C0, sorted_V, 'o', 'Color', current_color, ...
            'LineWidth', 1, 'MarkerSize', 5, 'DisplayName', sprintf('ΔT=%.1f K', target_deltaT));
        plot(sorted_C0, sorted_V, '-', 'Color', current_color, ...
            'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
end

set(gca, 'YScale', 'log');
xlabel('C₀ (wt%)');
ylabel('V (m s⁻¹)');
title(sprintf('%s - Velocity vs Concentration (Log Scale)', material_name));
legend('Location', 'best');
grid on;

% Continue with existing plots (Figure 6, 7, etc.)
create_additional_plots(all_C0_final_results, C0_values, concentration_colors, material_name);
create_3d_plots(all_C0_final_results, C0_values, material_name);
end



function create_additional_plots(all_C0_final_results, C0_values, concentration_colors, material_name)
%% CREATE_ADDITIONAL_PLOTS - Generate supplementary 2D visualizations
%
% PURPOSE:
% Creates additional 2D plots focusing on undercooling-dependent behaviour
% with concentration-differentiated visualization for detailed analysis.
%
% PLOT SPECIFICATIONS:
% Figure 6: V vs ΔT (Linear scale) - Shows velocity evolution with undercooling
% Figure 7: R vs ΔT (Log scale) - Emphasizes radius variation across wide range
%
% DATA PROCESSING:
% - Sorts data by undercooling for continuous line plotting
% - Uses filled markers with connecting lines for clarity
% - Applies consistent concentration-based color scheme
% - Handles missing data points gracefully

num_C0 = length(C0_values);

% 重新计算颜色映射 (与finalize_plots中保持一致)
if num_C0 > 1
    % 获取浓度的最小值和最大值
    C0_min = min(C0_values);
    C0_max = max(C0_values);

    % 为每个浓度计算其在范围内的相对位置 (0到1之间)
    C0_normalized = (C0_values - C0_min) / (C0_max - C0_min);

    % 获取完整的jet颜色映射表
    jet_colormap = jet(256);

    % 根据归一化的浓度位置插值得到对应颜色
    concentration_colors = zeros(num_C0, 3);
    for i = 1:num_C0
        % 将归一化位置映射到颜色表索引 (1到256)
        color_index = round(C0_normalized(i) * 255) + 1;
        color_index = max(1, min(256, color_index)); % 确保索引在有效范围内
        concentration_colors(i, :) = jet_colormap(color_index, :);
    end
else
    % 单一浓度用蓝色
    concentration_colors = [0, 0, 1];
end


% V with linear scale (Figure 6)
figure(6);
clf; hold on;
for c0_idx = 1:num_C0
    final_result = all_C0_final_results{c0_idx};
    if ~isempty(final_result) && size(final_result, 1) > 1
        current_color = concentration_colors(c0_idx, :);
        C0 = C0_values(c0_idx);

        % 排序数据以便绘制直线
        [sorted_deltaT, sort_idx] = sort(final_result(:,1));
        sorted_V = final_result(sort_idx, 2);

        % 绘制实心散点和直线连接
        plot(sorted_deltaT, sorted_V, '.', 'Color', current_color, ...
            'LineWidth', 2, 'MarkerSize', 15, 'DisplayName', sprintf('C0=%.1f', C0));
        plot(sorted_deltaT, sorted_V, '-', 'Color', current_color, ...
            'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
end
xlabel('ΔT (K)');
ylabel('V (m s^{-1})');
title(sprintf('%s - Velocity vs Undercooling (All C0 values) - Linear Scale', material_name));
legend('Location', 'best');
grid on;

% R with log scale (Figure 7)
figure(7);
clf; hold on;
for c0_idx = 1:num_C0
    final_result = all_C0_final_results{c0_idx};
    if ~isempty(final_result) && size(final_result, 1) > 1
        current_color = concentration_colors(c0_idx, :);
        C0 = C0_values(c0_idx);

        % 排序数据以便绘制直线
        [sorted_deltaT, sort_idx] = sort(final_result(:,1));
        sorted_R = final_result(sort_idx, 3);

        % 绘制实心散点和直线连接
        plot(sorted_deltaT, sorted_R, '.', 'Color', current_color, ...
            'LineWidth', 2, 'MarkerSize', 15, 'DisplayName', sprintf('C0=%.1f', C0));
        plot(sorted_deltaT, sorted_R, '-', 'Color', current_color, ...
            'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
end
% set(gca, 'YScale', 'log');
xlabel('ΔT (K)');
ylabel('R (m)');
title(sprintf('%s - Radius vs Undercooling (All C0 values) - Log Scale', material_name));
legend('Location', 'best');
grid on;
end




function create_3d_plots(all_C0_final_results, C0_values, material_name)
%% CREATE_3D_PLOTS - Generate three-dimensional visualization and surface fitting
%
% PURPOSE:
% Creates 3D scatter plots and performs comprehensive surface fitting
% analysis to reveal underlying relationships in the LGK model results.
%
% 3D VISUALIZATION FEATURES:
% Figure 4: V(C₀, ΔT) scatter plot with logarithmic Z-axis
% Figure 5: R(C₀, ΔT) scatter plot with logarithmic Z-axis
% - Color-coded by dependent variable magnitude
% - Optimized viewing angles for data clarity
% - Logarithmic scaling for wide dynamic ranges
%
% SURFACE FITTING WORKFLOW:
% 1. Data validation and outlier filtering
% 2. Multi-parameter regression analysis
% 3. Statistical quality assessment (R-squared values)
% 4. Custom functional form fitting for V and R
%
% MATHEMATICAL MODELS FITTED:
% - General form: ΔT = a×C₀^b×V^c (transport relationship)
% - V-specific: V = (ΔT/(a×C₀^b))^c (LGK velocity form)
% - R-specific: R = d×C₀^e×ΔT^f (LGK radius form)

fprintf('\nCreating 3D visualizations...\n');

% Prepare data for 3D plotting
all_C0_3D = [];
all_deltaT_3D = [];
all_V_3D = [];
all_R_3D = [];

num_C0 = length(C0_values);
for c0_idx = 1:num_C0
    C0 = C0_values(c0_idx);
    final_result = all_C0_final_results{c0_idx};

    if ~isempty(final_result)
        n_points = size(final_result, 1);
        all_C0_3D = [all_C0_3D; repmat(C0, n_points, 1)];
        all_deltaT_3D = [all_deltaT_3D; final_result(:, 1)];
        all_V_3D = [all_V_3D; final_result(:, 2)];
        all_R_3D = [all_R_3D; final_result(:, 3)];
    end
end

% Only create 3D plots if we have data
if ~isempty(all_C0_3D)
    % 3D scatter plot for V - Log Z scale
    figure(4);
    scatter3(all_C0_3D, all_deltaT_3D, all_V_3D, 20, all_V_3D, 'filled');
    xlabel('C0 (wt%)');
    ylabel('ΔT (K)');
    zlabel('V (m s^{-1})');
    title(sprintf('%s - 3D Visualization: Velocity (Log Z-axis)', material_name));
    set(gca, 'ZScale', 'log');
    colorbar;
    grid on;
    view(45, 30);

    % 3D scatter plot for R - Log Z scale
    figure(5);
    scatter3(all_C0_3D, all_deltaT_3D, all_R_3D, 20, all_R_3D, 'filled');
    xlabel('C0 (wt%)');
    ylabel('ΔT (K)');
    zlabel('R (m)');
    title(sprintf('%s - 3D Visualization: Radius (Log Z-axis)', material_name));
    set(gca, 'ZScale', 'log');
    colorbar;
    grid on;
    %     view(45, 30);
    view(135, 30);

    %     % 3D scatter plot for V - Linear scale (will be updated with surface fitting)
    %     % This will be created in the custom fitting function
    %     figure(8);
    %     scatter3(all_C0_3D, all_deltaT_3D, all_V_3D, 50, all_V_3D, 'filled');
    %     xlabel('C0 (wt%)');
    %     ylabel('ΔT (K)');
    %     zlabel('V (m s^{-1})');
    %     title('3D Visualization: Velocity (Linear Scale)');
    %     colorbar;
    %     grid on;
    %     view(45, 30);
    %
    %     % 3D scatter plot for R - Linear scale (will be updated with surface fitting)
    %     % This will be created in the custom fitting function
    %     figure(9);
    %     scatter3(all_C0_3D, all_deltaT_3D, all_R_3D, 50, all_R_3D, 'filled');
    %     xlabel('C0 (wt%)');
    %     ylabel('ΔT (K)');
    %     zlabel('R (m)');
    %     title('3D Visualization: Radius (Linear Scale)');
    %     colorbar;
    %     grid on;
    %     view(45, 30);

    % Surface fitting and summary statistics
    perform_surface_fitting(all_C0_3D, all_V_3D, all_deltaT_3D, all_R_3D, material_name);
    % Perform custom surface fitting for V and R with specified forms
    perform_custom_VR_fitting(all_C0_3D, all_deltaT_3D, all_V_3D, all_R_3D, material_name);
end

end

function perform_surface_fitting(all_C0_3D, all_V_3D, all_deltaT_3D, all_R_3D, material_name)
% Perform surface fitting and display statistics

% Summary statistics
fprintf('\n=== Summary Statistics ===\n');
total_points = length(all_V_3D);
fprintf('Total data points across all C0 values: %d\n', total_points);
if total_points > 0
    fprintf('C0 range: %.1f - %.1f wt%%\n', min(all_C0_3D), max(all_C0_3D));
    fprintf('ΔT range: %.2f - %.1f K\n', min(all_deltaT_3D), max(all_deltaT_3D));
    fprintf('V range: %.2e - %.2e m/s\n', min(all_V_3D), max(all_V_3D));
    fprintf('R range: %.2e - %.2e m\n', min(all_R_3D), max(all_R_3D));
else
    fprintf('No data points available for statistics\n');
    return;
end

% Surface fitting: ΔT = a * C0^b * V^c
fprintf('\n=== Surface Fitting: ΔT = a * C0^b * V^c ===\n');

if total_points >= 6  % Need at least 6 points for 3-parameter fitting
    % Prepare data for fitting
    C0_fit = all_C0_3D;
    V_fit = all_V_3D;
    deltaT_fit = all_deltaT_3D;

    % Remove any invalid data points
    valid_idx = ~isnan(C0_fit) & ~isnan(V_fit) & ~isnan(deltaT_fit) & ...
        C0_fit > 0 & V_fit > 0 & deltaT_fit > 0;
    C0_valid = C0_fit(valid_idx);
    V_valid = V_fit(valid_idx);
    deltaT_valid = deltaT_fit(valid_idx);

    if length(C0_valid) >= 6
        try
            % Logarithmic transformation: log(ΔT) = log(a) + b*log(C0) + c*log(V)
            log_C0 = log(C0_valid);
            log_V = log(V_valid);
            log_deltaT = log(deltaT_valid);

            % Create design matrix for multiple linear regression
            X = [ones(length(log_C0), 1), log_C0, log_V];

            % Solve using least squares: log_deltaT = X * coeffs
            coeffs = X \ log_deltaT;

            % Extract fitting parameters
            log_a = coeffs(1);
            b = coeffs(2);
            c = coeffs(3);
            a = exp(log_a);

            fprintf('Fitting successful!\n');
            fprintf('Fitted equation: ΔT = %.6e * C0^%.6f * V^%.6f\n', a, b, c);

            % Calculate R-squared
            deltaT_predicted = a * (C0_valid .^ b) .* (V_valid .^ c);
            SS_res = sum((deltaT_valid - deltaT_predicted).^2);
            SS_tot = sum((deltaT_valid - mean(deltaT_valid)).^2);
            R_squared = 1 - SS_res / SS_tot;
            fprintf('R-squared: %.6f\n', R_squared);

            % Create surface plot
            create_surface_plot(C0_valid, V_valid, deltaT_valid, a, b, c, R_squared, material_name);

        catch ME
            fprintf('Surface fitting failed: %s\n', ME.message);
            create_scatter_plot(C0_valid, V_valid, deltaT_valid, material_name);
        end
    else
        fprintf('Not enough valid data points for surface fitting (%d < 6)\n', length(C0_valid));
    end
else
    fprintf('Not enough total data points for surface fitting (%d < 6)\n', total_points);
end

end

function perform_custom_VR_fitting(all_C0_3D, all_deltaT_3D, all_V_3D, all_R_3D, material_name)
%% PERFORM_CUSTOM_VR_FITTING - Execute specialized V and R surface fitting
%
% PURPOSE:
% Performs surface fitting using the specific functional forms predicted
% by LGK theory, providing direct comparison with theoretical expectations.
%
% THEORETICAL BACKGROUND:
% LGK model predicts specific functional dependencies:
% - Velocity: V = (ΔT/(a×C₀^b))^c
% - Radius: R = d×C₀^e×ΔT^f
% These forms emerge from coupled diffusion and stability analysis.
%
% FITTING METHODOLOGY:
% V fitting: Logarithmic transformation to linear regression
% log(V) = c×log(ΔT) - c×b×log(C₀) - c×log(a)
% 
% R fitting: Direct logarithmic regression
% log(R) = log(d) + e×log(C₀) + f×log(ΔT)
%
% QUALITY METRICS:
% - R-squared correlation coefficients
% - Residual analysis for model validation
% - Parameter confidence assessment
%
% OUTPUTS:
% Fitted parameters stored in base workspace for Excel export:
% - V_surface_fitting: Contains V model parameters and statistics
% - R_surface_fitting: Contains R model parameters and statistics

if isempty(all_C0_3D)
    fprintf('No data available for custom V/R fitting\n');
    return;
end

fprintf('\n=== Custom V and R Surface Fitting ===\n');

% % Remove invalid data points
% valid_idx = ~isnan(all_C0_3D) & ~isnan(all_deltaT_3D) & ~isnan(all_V_3D) & ~isnan(all_R_3D) & ...
%     all_C0_3D > 0 & all_deltaT_3D > 0 & all_V_3D > 0 & all_R_3D > 0;

% Remove invalid data points with stricter filtering
valid_idx = ~isnan(all_C0_3D) & ~isnan(all_deltaT_3D) & ~isnan(all_V_3D) & ~isnan(all_R_3D) & ...
    all_C0_3D > 0 & all_deltaT_3D > 0 & all_V_3D > 0 & all_R_3D > 0 & ...
    isfinite(all_C0_3D) & isfinite(all_deltaT_3D) & isfinite(all_V_3D) & isfinite(all_R_3D) & ...
    all_R_3D < 1e-3 & all_R_3D > 1e-6;

C0_valid = all_C0_3D(valid_idx);
deltaT_valid = all_deltaT_3D(valid_idx);
V_valid = all_V_3D(valid_idx);
R_valid = all_R_3D(valid_idx);

if length(C0_valid) < 6
    fprintf('Not enough valid data points for custom fitting (%d < 6)\n', length(C0_valid));
    return;
end

% V fitting: V = (deltaT/(a*C0^b))^c
% Taking log: log(V) = c*log(deltaT) - c*log(a) - c*b*log(C0)
% Rearranging: log(V) = c*log(deltaT) - c*b*log(C0) - c*log(a)
% Linear form: log(V) = p1*log(deltaT) + p2*log(C0) + p3
% where p1 = c, p2 = -c*b, p3 = -c*log(a)
fprintf('\n--- V Fitting: V = (ΔT/(a×C0^b))^c ---\n');

try
    log_V = log(V_valid);
    log_deltaT = log(deltaT_valid);
    log_C0 = log(C0_valid);

    % Create design matrix
    X_V = [log_deltaT, log_C0, ones(length(log_V), 1)];

    % Solve linear system
    coeffs_V = X_V \ log_V;

    % Extract parameters
    c_V = coeffs_V(1);
    b_V = -coeffs_V(2) / c_V;
    a_V = exp(-coeffs_V(3) / c_V);

    % Calculate R-squared
    V_predicted = (deltaT_valid ./ (a_V * C0_valid.^b_V)).^c_V;
    SS_res_V = sum((V_valid - V_predicted).^2);
    SS_tot_V = sum((V_valid - mean(V_valid)).^2);
    R_squared_V = 1 - SS_res_V / SS_tot_V;

    fprintf('V fitting successful!\n');
    fprintf('V = (ΔT/(%.6e × C0^%.6f))^%.6f\n', a_V, b_V, c_V);
    fprintf('R-squared: %.6f\n', R_squared_V);

    % Create V surface plot
    create_V_surface_plot(C0_valid, deltaT_valid, V_valid, a_V, b_V, c_V, R_squared_V, material_name);

catch ME
    % Stop and cleanup progress timer if it exists
    if exist('progress_timer', 'var') && isvalid(progress_timer)
        stop(progress_timer);
        delete(progress_timer);
        fprintf('Progress timer stopped due to error.\n');
    end

    % Re-throw the error

    fprintf('V fitting failed: %s\n', ME.message);
end

% R fitting: R = d*C0^e*deltaT^f
% Taking log: log(R) = log(d) + e*log(C0) + f*log(deltaT)
% Linear form: log(R) = p1 + p2*log(C0) + p3*log(deltaT)
% where p1 = log(d), p2 = e, p3 = f
fprintf('\n--- R Fitting: R = d×C0^e×ΔT^f ---\n');

try
    log_R = log(R_valid);
    log_deltaT = log(deltaT_valid);
    log_C0 = log(C0_valid);

    % Create design matrix
    X_R = [ones(length(log_R), 1), log_C0, log_deltaT];

    % Solve linear system
    coeffs_R = X_R \ log_R;

    % Extract parameters
    d_R = exp(coeffs_R(1));
    e_R = coeffs_R(2);
    f_R = coeffs_R(3);

    % Calculate R-squared
    R_predicted = d_R * (C0_valid.^e_R) .* (deltaT_valid.^f_R);
    SS_res_R = sum((R_valid - R_predicted).^2);
    SS_tot_R = sum((R_valid - mean(R_valid)).^2);
    R_squared_R = 1 - SS_res_R / SS_tot_R;

    fprintf('R fitting successful!\n');
    fprintf('R = %.6e × C0^%.6f × ΔT^%.6f\n', d_R, e_R, f_R);
    fprintf('R-squared: %.6f\n', R_squared_R);

    % Create R surface plot
    create_R_surface_plot(C0_valid, deltaT_valid, R_valid, d_R, e_R, f_R, R_squared_R, material_name);

catch ME
    fprintf('R fitting failed: %s\n', ME.message);
end

end

function create_V_surface_plot(C0_valid, deltaT_valid, V_valid, a_V, b_V, c_V, R_squared_V, material_name)
% Create V surface plot with fitted equation: V = (deltaT/(a*C0^b))^c

figure(8);
clf;

% Plot original data points
scatter3(C0_valid, deltaT_valid, V_valid, 10, V_valid, 'filled', 'DisplayName', 'Data Points');
hold on;

% Create fitted surface
C0_grid = linspace(min(C0_valid), max(C0_valid), 30);
deltaT_grid = linspace(min(deltaT_valid), max(deltaT_valid), 30);
[C0_mesh, deltaT_mesh] = meshgrid(C0_grid, deltaT_grid);
V_mesh = (deltaT_mesh ./ (a_V * C0_mesh.^b_V)).^c_V;

% Plot fitted surface with proper color mapping based on Z-axis (V values)
surf(C0_mesh, deltaT_mesh, V_mesh, V_mesh, 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Fitted Surface');

% Formatting
xlabel('C0 (wt%)');
ylabel('ΔT (K)');
zlabel('V (m s^{-1})');
title(sprintf('%s - V Fit: V = (ΔT/(%.6g×C0^{%.6g}))^{%.6g} (R² = %.4f)', material_name, a_V, b_V, c_V, R_squared_V));
colormap(parula);
shading interp;  % Smooth gradient shading to eliminate grid effects
caxis([min(V_mesh(:)), max(V_mesh(:))]);  % Set the colour range based on Z-axis values
colorbar;
grid on;
view(45, 30);
legend('Location', 'best');

% Add text annotation with equation in top-left corner using annotation
text_str = sprintf('V = (ΔT/(%.6g×C0^{%.6g}))^{%.6g}\nR² = %.4f', a_V, b_V, c_V, R_squared_V);
annotation('textbox', [0.05, 0.81, 0.3, 0.1], 'String', text_str, ...
    'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'FitBoxToText', 'on');

% Save V fitting results to base workspace for Excel export
V_surface_fitting = struct();
V_surface_fitting.equation = 'V = (ΔT/(a×C0^b))^c';
V_surface_fitting.a = a_V;
V_surface_fitting.b = b_V;
V_surface_fitting.c = c_V;
V_surface_fitting.R_squared = R_squared_V;
V_surface_fitting.figure_number = 8;
assignin('base', 'V_surface_fitting', V_surface_fitting);

end

function create_R_surface_plot(C0_valid, deltaT_valid, R_valid, d_R, e_R, f_R, R_squared_R, material_name)
% Create R surface plot with fitted equation: R = d*C0^e*deltaT^f

figure(9);
clf;

% Plot original data points
scatter3(C0_valid, deltaT_valid, R_valid, 10, R_valid, 'filled', 'DisplayName', 'Data Points');
hold on;

% Create fitted surface
C0_grid = linspace(min(C0_valid), max(C0_valid), 30);
deltaT_grid = linspace(min(deltaT_valid), max(deltaT_valid), 30);
[C0_mesh, deltaT_mesh] = meshgrid(C0_grid, deltaT_grid);
R_mesh = d_R * (C0_mesh.^e_R) .* (deltaT_mesh.^f_R);

% Plot fitted surface with proper color mapping based on Z-axis (R values)
surf(C0_mesh, deltaT_mesh, R_mesh, R_mesh, 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Fitted Surface');

% Formatting
xlabel('C0 (wt%)');
ylabel('ΔT (K)');
zlabel('R (m)');
title(sprintf('%s - R Fit: R = %.6g×C0^{%.6g}×ΔT^{%.6g} (R² = %.4f)', material_name, d_R, e_R, f_R, R_squared_R));
colormap(parula);
shading interp;  % Smooth gradient shading to eliminate grid effects
caxis([min(R_mesh(:)), max(R_mesh(:))]);  % Set the colour range based on Z-axis values
colorbar;
grid on;
view(135, 30);
legend('Location', 'best');

% Add text annotation with equation in top-left corner using annotation
text_str = sprintf('R = %.6g×C0^{%.6g}×ΔT^{%.6g}\nR² = %.4f', d_R, e_R, f_R, R_squared_R);
annotation('textbox', [0.05, 0.81, 0.3, 0.1], 'String', text_str, ...
    'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'FitBoxToText', 'on');

% Save R fitting results to base workspace for Excel export
R_surface_fitting = struct();
R_surface_fitting.equation = 'R = d×C0^e×ΔT^f';
R_surface_fitting.d = d_R;
R_surface_fitting.e = e_R;
R_surface_fitting.f = f_R;
R_surface_fitting.R_squared = R_squared_R;
R_surface_fitting.figure_number = 9;
assignin('base', 'R_surface_fitting', R_surface_fitting);

end

function create_surface_plot(C0_valid, V_valid, deltaT_valid, a, b, c, R_squared, material_name)
% Create surface plot with fitted equation

figure(10);

% Plot original data points
scatter3(C0_valid, V_valid, deltaT_valid, 10, deltaT_valid, 'filled', 'DisplayName', 'Data Points');
hold on;

% Create fitted surface
C0_grid = linspace(min(C0_valid), max(C0_valid), 30);
V_grid = logspace(log10(min(V_valid)), log10(max(V_valid)), 30);
[C0_mesh, V_mesh] = meshgrid(C0_grid, V_grid);
deltaT_mesh = a * (C0_mesh .^ b) .* (V_mesh .^ c);

% Plot fitted surface with proper color mapping based on Z-axis (deltaT values)
surf(C0_mesh, V_mesh, deltaT_mesh, deltaT_mesh, 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Fitted Surface');

% Formatting
xlabel('C0 (wt%)');
ylabel('V (m s^{-1})');
zlabel('ΔT (K)');
title(sprintf('%s - Surface Fit: ΔT = %.6g × C0^{%.6g} × V^{%.6g} (R² = %.4f)', material_name, a, b, c, R_squared));
% set(gca, 'YScale', 'log');
colormap(parula);
shading interp;  % Smooth gradient shading to eliminate grid effects
caxis([min(deltaT_mesh(:)), max(deltaT_mesh(:))]);  % Set the colour range based on Z-axis values
colorbar;
grid on;
view(45, 30);
legend('Location', 'best');

% Add text annotation with equation in top-left corner using annotation
text_str = sprintf('ΔT = %.6g × C0^{%.6g} × V^{%.6g}\nR² = %.4f', a, b, c, R_squared);
annotation('textbox', [0.05, 0.81, 0.3, 0.1], 'String', text_str, ...
    'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'FitBoxToText', 'on');

% Save surface fitting results to base workspace for Excel export
deltaT_surface_fitting = struct();
deltaT_surface_fitting.equation = 'ΔT = a×C0^b×V^c';
deltaT_surface_fitting.a = a;
deltaT_surface_fitting.b = b;
deltaT_surface_fitting.c = c;
deltaT_surface_fitting.R_squared = R_squared;
deltaT_surface_fitting.figure_number = 10;
assignin('base', 'deltaT_surface_fitting', deltaT_surface_fitting);

end

function create_scatter_plot(C0_valid, V_valid, deltaT_valid, material_name)
% Create scatter plot when fitting fails

figure(10);
scatter3(C0_valid, V_valid, deltaT_valid, 10, deltaT_valid, 'filled');
xlabel('C0 (wt%)');
ylabel('V (m s^{-1})');
zlabel('ΔT (K)');
title(sprintf('%s - Data Points (Surface Fitting Failed)', material_name));
set(gca, 'YScale', 'log');
colorbar;
grid on;
view(45, 30);

end


function adjust_figure_layout()
%% ADJUST_FIGURE_LAYOUT - Optimize figure positioning for comprehensive display
%
% PURPOSE:
% Implements intelligent figure layout management for optimal screen
% utilization and scientific presentation quality.
%
% LAYOUT STRATEGY:
% Two-row arrangement optimized for widescreen displays:
% Row 1 (top): Figures 1, 2, 6, 7 - 2D analysis plots
% Row 2 (bottom): Figures 4, 5, 8, 9, 10 - 3D and surface fitting plots
%
% POSITIONING ALGORITHM:
% - Automatic screen size detection and adaptation
% - Uniform horizontal spacing with aspect ratio preservation
% - Vertical positioning for non-overlapping display
% - Figure size optimization for readability and detail preservation
%
% SCIENTIFIC PRESENTATION FEATURES:
% - Consistent figure dimensions for publication quality
% - Optimal viewing angles for 3D plots
% - Coordinated color schemes across related visualizations
% - Professional layout suitable for technical documentation

% Get screen size
screenSize = get(0, 'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4);

% Calculate figure dimensions (half screen height, keep aspect ratio)
figHeight = screenHeight * 0.5;
figWidth = figHeight * 1.3;  % Maintain aspect ratio

% Calculate positions
% First row: Figure 1, 2, 6, 7 - bottom at 40% of screen height
% Second row: Figure 4, 5, 8, 9, 10 - bottom at 10% of screen height
firstRowY = screenHeight * 0.4;   % Bottom at 40% of screen height
secondRowY = screenHeight * 0;  % Bottom at 10% of screen height

% Calculate horizontal positions for uniform distribution
% First row: 4 figures (1, 2, 6, 7)
firstRowFigures = [1, 2, 6, 7];
firstRowCount = length(firstRowFigures);
firstRowTotalWidth = firstRowCount * figWidth;
firstRowSpacing = (screenWidth - figWidth) / (firstRowCount - 1);

% Second row: 5 figures (4, 5, 8, 9, 10)
secondRowFigures = [4, 5, 8, 9, 10];
secondRowCount = length(secondRowFigures);
secondRowTotalWidth = secondRowCount * figWidth;
secondRowSpacing = (screenWidth - figWidth) / (secondRowCount - 1);


fprintf('Figure layout parameters:\n');
fprintf('Screen size: %.0f x %.0f\n', screenWidth, screenHeight);
fprintf('Figure size: %.0f x %.0f\n', figWidth, figHeight);
fprintf('First row Y position: %.0f\n', firstRowY);
fprintf('Second row Y position: %.0f\n', secondRowY);
fprintf('First row spacing: %.0f\n', firstRowSpacing);
fprintf('Second row spacing: %.0f\n', secondRowSpacing);

% Position first row figures (1, 2, 6, 7)
for i = 1:length(firstRowFigures)
    figNum = firstRowFigures(i);

    if ishandle(figNum)
        figure(figNum);

        % Calculate x position for uniform distribution
        % Figure 1 starts at x=0, Figure 7 ends at x=screenWidth
        xPos = (i - 1) * firstRowSpacing;

        % Set figure position
        set(gcf, 'Position', [xPos, firstRowY, figWidth, figHeight]);

        fprintf('First row - Figure %d position: [%.0f, %.0f, %.0f, %.0f]\n', ...
            figNum, xPos, firstRowY, figWidth, figHeight);
    else
        fprintf('Warning: Figure %d does not exist, skipping.\n', figNum);
    end
end

% Position second row figures (4, 5, 8, 9, 10)
for i = 1:length(secondRowFigures)
    figNum = secondRowFigures(i);

    if ishandle(figNum)
        figure(figNum);

        % Calculate x position for uniform distribution
        % Figure 4 starts at x=0, Figure 10 ends at x=screenWidth
        xPos = (i - 1) * secondRowSpacing;

        % Set figure position
        set(gcf, 'Position', [xPos, secondRowY, figWidth, figHeight]);

        fprintf('Second row - Figure %d position: [%.0f, %.0f, %.0f, %.0f]\n', ...
            figNum, xPos, secondRowY, figWidth, figHeight);
    else
        fprintf('Warning: Figure %d does not exist, skipping.\n', figNum);
    end
end

% Bring control panel (Figure 999) to top if it exists
if ishandle(999)
    figure(999);
    fprintf('Control panel brought to top\n');
end

fprintf('Figure layout adjustment completed.\n');
end





function V_guess = get_SCN_fitted_V_guess(C0, deltaT)
% V_guess function based on sixth-degree polynomial fitting

V_05K_poly_coeffs = [-6.659290507617004e-04, 2.369698931357119e-03, -3.396311141455152e-03, 2.487586100839964e-03, -9.373914269260196e-04, 1.173158416800813e-04, 3.764308570717425e-05];
V_09K_poly_coeffs = [-5.533479076924480e-03, 1.958381856164818e-02, -2.749709894968078e-02, 1.941241773272147e-02, -7.029175352673138e-03, 9.566442665513241e-04, 1.691807667749731e-04];

V_at_05K = polyval(V_05K_poly_coeffs, C0);
V_at_09K = polyval(V_09K_poly_coeffs, C0);
deltaT_clamped = max(0.5, min(0.9, deltaT));
weight_05K = (0.9 - deltaT_clamped) / 0.4;
weight_09K = (deltaT_clamped - 0.5) / 0.4;
V_guess = weight_05K * V_at_05K + weight_09K * V_at_09K;
V_guess = max(V_guess, 1e-10);
end

function R_guess = get_SCN_fitted_R_guess(C0, deltaT)
% R_guess function based on sixth-degree polynomial fitting

R_05K_poly_coeffs = [3.933800623029993e-04, -1.399829155419792e-03, 1.993583146697738e-03, -1.455094212818860e-03, 5.807424879235848e-04, -1.226741862485822e-04, 2.463444587967703e-05];
R_09K_poly_coeffs = [4.062364954450699e-04, -1.402966115638020e-03, 1.906008843886529e-03, -1.294100697718033e-03, 4.630895797590333e-04, -8.404160371893581e-05, 1.193829099220079e-05];

R_at_05K = polyval(R_05K_poly_coeffs, C0);
R_at_09K = polyval(R_09K_poly_coeffs, C0);
deltaT_clamped = max(0.5, min(0.9, deltaT));
weight_05K = (0.9 - deltaT_clamped) / 0.4;
weight_09K = (deltaT_clamped - 0.5) / 0.4;
R_guess = weight_05K * R_at_05K + weight_09K * R_at_09K;
R_guess = max(R_guess, 1e-10);
end



% PROCESS_SINGLE_C0_WORKER - Process single C0 value in parallel worker
function [worker_results, worker_final, worker_min_dt] = process_single_c0_worker(...
    c0_idx, C0, material_params, deltaT_sequence, max_iterations, ...
    V_a, V_b, V_c, R_d, R_e, R_f, use_analytical_approximation)
%% PROCESS_SINGLE_C0_WORKER - Parallel worker function for single concentration
%
% PURPOSE:
% Executes complete LGK calculation workflow for a single concentration
% value within parallel computing environment, handling material-specific
% computational strategies and parameter optimization.
%
% WORKER FUNCTIONALITY:
% 1. Material type detection (SCN vs standard alloys)
% 2. Initial guess function configuration (analytical vs parametric)
% 3. Iterative calculation with parameter refinement
% 4. Convergence monitoring and result validation
%
% MATERIAL-SPECIFIC STRATEGIES:
% SCN-Acetone: Uses polynomial interpolation for specialized deltaT range
% Standard alloys: Employs parametric or analytical initial guess functions
% Analytical mode: Utilizes theoretical equations (8.91, 8.92) when available
%
% PARALLEL OPTIMIZATION:
% - Minimal inter-worker communication requirements
% - Self-contained error handling and recovery
% - Efficient memory management for large parameter spaces
% - Progress reporting without synchronization overhead
%
% CONVERGENCE STRATEGY:
% - Multi-scale initial guess testing with scaling factors
% - Adaptive parameter refinement based on convergence history
% - Early termination for non-convergent cases
% - Quality-based result selection and inheritance

% Debug: Check if material_params are correctly passed to worker
if use_analytical_approximation
    fprintf('Worker %d: material_params check for C0=%.6f:\n', c0_idx, C0);
    fprintf('  Gibbs_Tom = %.6e\n', material_params.Gibbs_Tom);
    fprintf('  m_Al = %.6f\n', material_params.m_Al);
    fprintf('  k0_Al = %.6f\n', material_params.k0_Al);
    fprintf('  D_Al = %.6e\n', material_params.D_Al);
    
    % Check if analytical function exists in worker
    if exist('create_analytical_guess_functions', 'file') ~= 2
        fprintf('ERROR: create_analytical_guess_functions not found in worker %d\n', c0_idx);
        % Try to add path again
        current_path = pwd;
        addpath(current_path);
        fprintf('Added path %s to worker %d\n', current_path, c0_idx);
    end
end

% Check if this is SCN material based on deltaT sequence
is_SCN_material = (length(deltaT_sequence) == 2 && deltaT_sequence(1) == 0.5 && deltaT_sequence(2) == 0.9);

if is_SCN_material
    % For SCN material, use polynomial interpolation
    fprintf('Worker %d: Using SCN polynomial interpolation for C0 = %.1f\n', c0_idx, C0);
    
    % Create guess functions using polynomial interpolation
    V_guess_func = @(deltaT) arrayfun(@(dt) get_SCN_fitted_V_guess(C0, dt), deltaT);
    R_guess_func = @(deltaT) arrayfun(@(dt) get_SCN_fitted_R_guess(C0, dt), deltaT);
    
    % Calculate directly for the two deltaT values
    [result, converged, last_diverged_deltaT] = calculate_VR(V_guess_func, R_guess_func, C0, deltaT_sequence, material_params);
    
    % Store results
    worker_results = {result};
    worker_final = result;
    if ~isempty(result)
        worker_min_dt = result(end, 1);
    else
        worker_min_dt = NaN;
    end
    
    fprintf('Worker %d: SCN calculation completed for C0 = %.1f with %d points\n', c0_idx, C0, size(result, 1));
    
else
    % For non-SCN materials, check if using analytical approximation
    % Check if analytical function is available
    if exist('create_analytical_guess_functions', 'file') ~= 2
        fprintf('Worker %d: Analytical function not available, falling back to iterative\n', c0_idx);
        use_analytical_approximation = false;
    end
    
    if use_analytical_approximation
        try
            % Create analytical guess functions with parameter validation
            fprintf('Worker %d: Calling create_analytical_guess_functions with:\n', c0_idx);
            fprintf('  C0=%.6f, Gibbs_Tom=%.6e, m_Al=%.6f, k0_Al=%.6f, D_Al=%.6e\n', ...
                    C0, material_params.Gibbs_Tom, material_params.m_Al, material_params.k0_Al, material_params.D_Al);
            
            [V_guess_func, R_guess_func] = create_analytical_guess_functions(C0, ...
                material_params.Gibbs_Tom, material_params.m_Al, material_params.k0_Al, material_params.D_Al);
            
            % Test the functions
            test_deltaT = 1.0;
            test_V = V_guess_func(test_deltaT);
            test_R = R_guess_func(test_deltaT);
            fprintf('Worker %d: Function test - V(%.1f)=%.6e, R(%.1f)=%.6e\n', ...
                    c0_idx, test_deltaT, test_V, test_deltaT, test_R);
            
        catch ME
            fprintf('ERROR in worker %d calling analytical functions: %s\n', c0_idx, ME.message);
            fprintf('Stack trace:\n');
            for i = 1:length(ME.stack)
                fprintf('  %s line %d\n', ME.stack(i).name, ME.stack(i).line);
            end
            
            % Fall back to error values
            worker_results = {[]};
            worker_final = [];
            worker_min_dt = NaN;
            return;
        end
        
        % Calculate directly using analytical functions
        [result, converged, last_diverged_deltaT] = calculate_VR(V_guess_func, R_guess_func, C0, deltaT_sequence, material_params);
        
        % Store results
        worker_results = {result};
        worker_final = result;
        if ~isempty(result)
            worker_min_dt = result(end, 1);
        else
            worker_min_dt = NaN;
        end
        
        fprintf('Worker %d: Analytical calculation completed for C0 = %.6f with %d points\n', c0_idx, C0, size(result, 1));
        return;
    else
        fprintf('Worker %d: Using iterative approach for C0 = %.1f\n', c0_idx, C0);
        
        % Initialize parameters
        V_a_current = V_a;
        V_b_current = V_b;
        V_c_current = V_c;
        R_d_current = R_d;
        R_e_current = R_e;
        R_f_current = R_f;
    end
    
    % Store iteration results
    all_iterations = cell(max_iterations, 1);
    min_cnvg_DT = zeros(max_iterations, 1);
    
    % Try different V_a scaling factors if needed
    scaling_factors = [1, 2, 0.5, 4, 0.25];
    best_result = [];
    
    for scale_idx = 1:length(scaling_factors)
        scale_factor = scaling_factors(scale_idx);
        V_a_test = V_a * scale_factor;
        
        % Create guess functions
        V_guess_func = @(deltaT) (deltaT./(V_a_test * C0.^V_b_current)).^V_c_current;
        R_guess_func = @(deltaT) R_d_current * (C0.^R_e_current) .* (deltaT.^R_f_current);
        
        % Test with first few deltaT values
        test_deltaT = deltaT_sequence(1:min(3, length(deltaT_sequence)));
        [test_result, ~, ~] = calculate_VR(V_guess_func, R_guess_func, C0, test_deltaT, material_params);
        
        if size(test_result, 1) >= 2
            best_result = test_result;
            V_a_current = V_a_test;
            break;
        end
    end
    
    if isempty(best_result)
        % No successful scaling found
        worker_results = {[]};
        worker_final = [];
        worker_min_dt = NaN;
        fprintf('Worker %d: No successful scaling found for C0 = %.1f\n', c0_idx, C0);
        return;
    end
    
    % Perform iterative calculation
    for iter = 1:max_iterations
        % Create guess functions
        V_guess_func = @(deltaT) (deltaT./(V_a_current * C0.^V_b_current)).^V_c_current;
        R_guess_func = @(deltaT) R_d_current * (C0.^R_e_current) .* (deltaT.^R_f_current);
        
        % Calculate
        [result, converged, last_diverged_deltaT] = calculate_VR(V_guess_func, R_guess_func, C0, deltaT_sequence, material_params);
        
        % Store results
        all_iterations{iter} = result;
        if ~isempty(result)
            min_cnvg_DT(iter) = result(end, 1);
        else
            min_cnvg_DT(iter) = NaN;
        end
        
        % Check convergence
        if size(result, 1) == length(deltaT_sequence) || size(result, 1) < 3
            break;
        end
        
        % Update parameters based on fitting (simplified for parallel worker)
        if size(result, 1) >= 3
            % Simple parameter update
            deltaT_data = result(:, 1);
            V_data = result(:, 2);
            R_data = result(:, 3);
            
            % Update V parameters
            log_V = log(V_data);
            log_deltaT = log(deltaT_data);
            log_C0 = log(repmat(C0, length(deltaT_data), 1));
            
            X_V = [log_deltaT, log_C0, ones(length(log_V), 1)];
            coeffs_V = X_V \ log_V;
            
            V_c_current = coeffs_V(1);
            V_b_current = -coeffs_V(2) / V_c_current;
            V_a_current = exp(-coeffs_V(3) / V_c_current);
            
            % Update R parameters
            log_R = log(R_data);
            X_R = [ones(length(log_R), 1), log_C0, log_deltaT];
            coeffs_R = X_R \ log_R;
            
            R_d_current = exp(coeffs_R(1));
            R_e_current = coeffs_R(2);
            R_f_current = coeffs_R(3);
        end
    end
    
    % Return results
    worker_results = all_iterations;
    worker_final = all_iterations{iter};
    worker_min_dt = min_cnvg_DT(1:iter);
    
    fprintf('Worker %d: Iterative calculation completed for C0 = %.1f with %d points\n', c0_idx, C0, size(worker_final, 1));
end

end


function progress_manager = init_progress_manager(material_name, total_points, gui_handles, pool_already_ready)
%% INIT_PROGRESS_MANAGER - Initialize comprehensive progress tracking system
%
% PURPOSE:
% Establishes sophisticated progress monitoring with historical data analysis,
% time estimation algorithms, and GUI synchronization for optimal user experience.
%
% PROGRESS TRACKING FEATURES:
% - Historical computation time analysis for accurate estimation
% - Parallel pool initialization time accounting
% - Real-time GUI updates with progress bars and time displays
% - Material-specific performance profiling
%
% TIME ESTIMATION ALGORITHM:
% Uses historical data when available:
% time_per_point = (total_time - parallel_init_time) / total_points
% estimated_total = parallel_init_time + (current_points × time_per_point)
%
% INPUTS:
%   material_name - Current material identifier for historical lookup
%   total_points - Expected total calculation points for progress scaling
%   gui_handles - GUI component handles for real-time updates
%   pool_already_ready - Boolean indicating parallel pool pre-initialization
%
% OUTPUTS:
%   progress_manager - Structure containing timing data and GUI references


% Initialise the progress manager
if nargin < 4
    pool_already_ready = false;
end

progress_manager = struct();
progress_manager.material_name = material_name;
progress_manager.total_points = total_points;
progress_manager.pool_already_ready = pool_already_ready;

% Validate input parameters
if isempty(material_name) || ~ischar(material_name)
    material_name = 'Unknown';
    progress_manager.material_name = material_name;
    fprintf('Warning: Invalid material name, using default\n');
end

if isempty(total_points) || ~isnumeric(total_points) || total_points <= 0
    total_points = 100;  % Default value
    progress_manager.total_points = total_points;
    fprintf('Warning: Invalid total_points, using default value: %d\n', total_points);
end

progress_manager.gui_handles = gui_handles;
progress_manager.start_time = now;
progress_manager.parallel_init_time = 0;
progress_manager.completed_points = 0;

% Load historical data
hist_data = load_timing_history(material_name);
progress_manager.hist_data = hist_data;

% Estimate time – adjust based on parallel pool status
if ~isempty(hist_data)
    % Estimate based on historical data
    time_per_point = (hist_data.total_time - hist_data.parallel_init_time) / hist_data.total_points;
    progress_manager.estimated_calc_time = time_per_point * total_points;
    
    % If the parallel pool is already ready, no initialisation time is needed
    if pool_already_ready
        progress_manager.estimated_parallel_time = 0;
        fprintf('Using historical data - pool ready, no init time needed\n');
    else
        progress_manager.estimated_parallel_time = min(hist_data.parallel_init_time, 0.5); % Maximum wait time of 0.5 minutes
        fprintf('Using historical data - pool may need wait time: %.2f min\n', progress_manager.estimated_parallel_time);
    end
    
    progress_manager.estimated_total_time = progress_manager.estimated_parallel_time + progress_manager.estimated_calc_time;
else
    % Default estimate when no historical data is available
    if pool_already_ready
        progress_manager.estimated_parallel_time = 0; % Pool is ready
        progress_manager.estimated_calc_time = total_points * 0.5 / 60; % 0.5s per point
        fprintf('No historical data - pool ready, estimating calc time only\n');
    else
        progress_manager.estimated_parallel_time = 0.5; % 0.5 seconds per point
        progress_manager.estimated_calc_time = total_points * 0.6 / 60; % Slightly increased time per point
        fprintf('No historical data - pool may need setup, conservative estimate\n');
    end
    progress_manager.estimated_total_time = progress_manager.estimated_parallel_time + progress_manager.estimated_calc_time;
end

fprintf('Progress Manager Initialized:\n');
fprintf('  Material: %s\n', material_name);
fprintf('  Total points: %d\n', total_points);
fprintf('  Pool already ready: %s\n', mat2str(pool_already_ready));
fprintf('  Estimated parallel init: %.2f min\n', progress_manager.estimated_parallel_time);
fprintf('  Estimated calculation: %.2f min\n', progress_manager.estimated_calc_time);
fprintf('  Estimated total: %.2f min\n', progress_manager.estimated_total_time);

% Initial display
update_progress_display(progress_manager, 0, 'Calculation Starting');
end

function hist_data = load_timing_history(material_name)
% Load historical timing data
hist_data = [];
filename = 'calculation_timing_history.mat';

% Clean the material name to make it a valid struct field name
clean_material_name = clean_field_name(material_name);

try
    if exist(filename, 'file')
        data = load(filename, 'timing_history');
        if isfield(data, 'timing_history') && isfield(data.timing_history, clean_material_name)
            hist_data = data.timing_history.(clean_material_name);
            fprintf('Loaded timing history for %s: %d points, %.2f min total\n', ...
                    material_name, hist_data.total_points, hist_data.total_time);
        else
            fprintf('No timing history found for material: %s\n', material_name);
        end
    else
        fprintf('No timing history file found\n');
    end
catch ME
    fprintf('Error loading timing history: %s\n', ME.message);
end
end

function save_timing_history(material_name, total_points, parallel_init_time, total_time)
% Save historical timing data
filename = 'calculation_timing_history.mat';

% Clean the material name to make it a valid struct field name
clean_material_name = clean_field_name(material_name);

try
    % Load existing data
    if exist(filename, 'file')
        data = load(filename, 'timing_history');
        timing_history = data.timing_history;
    else
        timing_history = struct();
    end
    
    % Update data for the current material
    timing_history.(clean_material_name).total_points = total_points;
    timing_history.(clean_material_name).parallel_init_time = parallel_init_time;
    timing_history.(clean_material_name).total_time = total_time;
    timing_history.(clean_material_name).last_updated = datestr(now);
    timing_history.(clean_material_name).original_name = material_name; % Save the original name
    timing_history.(clean_material_name).pool_was_ready = parallel_init_time < 0.1; % Record whether the pool was already ready
    
    % save
    save(filename, 'timing_history');
    fprintf('Timing history saved for %s: %d points, %.2f min parallel, %.2f min total\n', ...
            material_name, total_points, parallel_init_time, total_time);
    
catch ME
    fprintf('Error saving timing history: %s\n', ME.message);
end
end

function progress_manager = update_progress_checkpoint(progress_manager, checkpoint_type, completed_points, message)
%% UPDATE_PROGRESS_CHECKPOINT - Process major calculation milestones
%
% PURPOSE:
% Updates progress tracking at key computational milestones with intelligent
% time estimation and GUI feedback synchronization.
%
% CHECKPOINT TYPES:
% 'parallel_init_complete': Parallel pool initialization finished
% 'calculation_progress': Incremental calculation progress updates
% 'calculation_complete': Final completion with historical data storage
%
% PROGRESS CALCULATION STRATEGY:
% - Parallel initialization: 0-15% depending on pool status
% - Calculation progress: 15-99% based on completed points
% - Final completion: 100% with historical data archival
%
% GUI SYNCHRONIZATION:
% Updates progress bar, percentage display, elapsed time, and
% remaining time estimation with rate-limited refresh for performance

current_time = now;
elapsed_time = (current_time - progress_manager.start_time) * 24 * 60; % min

switch checkpoint_type
    case 'parallel_init_complete'
        progress_manager.parallel_init_time = elapsed_time;
        if progress_manager.pool_already_ready
            % If the pool was already ready, this time should be very short, with minimal progress
            progress_percent = min((elapsed_time / progress_manager.estimated_total_time) * 100, 3);
        else
            % If the pool needs to wait/initialise, the progress is slightly larger
            progress_percent = min((elapsed_time / progress_manager.estimated_total_time) * 100, 15);
        end
        
    case 'calculation_progress'
        progress_manager.completed_points = completed_points;
        % Calculate progress: parallel initialisation progress + calculation progress
        if isfield(progress_manager, 'parallel_init_time')
            parallel_progress = (progress_manager.parallel_init_time / progress_manager.estimated_total_time) * 100;
        else
            parallel_progress = 0;
        end
        calc_progress = (completed_points / progress_manager.total_points) * 85; % Calculation portion accounts for 85%
        progress_percent = parallel_progress + calc_progress;
            
    case 'calculation_complete'
        progress_manager.completed_points = progress_manager.total_points;
        progress_percent = 100;
        % Save historical data
        actual_parallel_time = progress_manager.parallel_init_time;
        if progress_manager.pool_already_ready && actual_parallel_time < 0.1
            actual_parallel_time = 0; % Pool already ready, record as 0 initialisation time
        end
        save_timing_history(progress_manager.material_name, progress_manager.total_points, ...
                           actual_parallel_time, elapsed_time);
        
    otherwise
        progress_percent = (elapsed_time / progress_manager.estimated_total_time) * 100;
end

% Limit progress range
progress_percent = max(0, min(99.9, progress_percent));
if strcmp(checkpoint_type, 'calculation_complete')
    progress_percent = 100;
end

% Update display
update_progress_display(progress_manager, progress_percent, message);

fprintf('Progress Checkpoint [%s]: %.1f%% | Elapsed: %.2f min | Message: %s\n', ...
        checkpoint_type, progress_percent, elapsed_time, message);
end

function update_progress_display(progress_manager, progress_percent, message)
% Update GUI progress display
try
    gui_handles = progress_manager.gui_handles;
    if isempty(gui_handles) || ~isstruct(gui_handles)
        return;
    end
    
    % Check whether GUI is still valid
    if isfield(gui_handles, 'figure') && ~isempty(gui_handles.figure)
        if ~isvalid(gui_handles.figure) || ~ishandle(gui_handles.figure)
            return;
        end
    end
    
    current_time = now;
    elapsed_time = (current_time - progress_manager.start_time) * 24 * 60; % minutes
    
    % Update progress bar
    if isfield(gui_handles, 'progress_text_bar') && isvalid(gui_handles.progress_text_bar)
        bar_length = 56; % Fixed length
        filled_chars = round(bar_length * progress_percent / 100);
        empty_chars = bar_length - filled_chars;
        progress_bar_str = ['[' repmat('█', 1, filled_chars) repmat('░', 1, empty_chars) ']'];
        set(gui_handles.progress_text_bar, 'String', progress_bar_str);
    end
    
    % Update percentage
    if isfield(gui_handles, 'progress_percent') && isvalid(gui_handles.progress_percent)
        set(gui_handles.progress_percent, 'String', sprintf('%.1f%%', progress_percent));
    end
    
    % Update elapsed time
    if isfield(gui_handles, 'progress_elapsed') && isvalid(gui_handles.progress_elapsed)
        if elapsed_time < 1
            elapsed_str = sprintf('Elapsed: %.0f sec', elapsed_time * 60);
        else
            elapsed_str = sprintf('Elapsed: %.2f min', elapsed_time);
        end
        set(gui_handles.progress_elapsed, 'String', elapsed_str);
    end
    
    % Update remaining time
    if isfield(gui_handles, 'progress_remaining') && isvalid(gui_handles.progress_remaining)
        if progress_percent >= 100
            remaining_str = 'Remaining: Complete';
        elseif progress_percent <= 0
            remaining_str = sprintf('Remaining: %.1f min', progress_manager.estimated_total_time);
        else
            % Display at least 0.1 minutes
            estimated_remaining = (elapsed_time / progress_percent * 100) - elapsed_time;
            estimated_remaining = max(0.1, estimated_remaining);
            
            if estimated_remaining < 1
                remaining_str = sprintf('Remaining: %.0f sec', estimated_remaining * 60);
            else
                remaining_str = sprintf('Remaining: %.1f min', estimated_remaining);
            end
        end
        set(gui_handles.progress_remaining, 'String', remaining_str);
    end
    
    % Force GUI update
    drawnow limitrate;
    
catch ME
    fprintf('Error updating progress display: %s\n', ME.message);
end
end

function clean_name = clean_field_name(material_name)
% Clean material name to make it a valid MATLAB structure field name
if isempty(material_name) || ~ischar(material_name)
    clean_name = 'Unknown_Material';
    return;
end

clean_name = material_name;
clean_name = strrep(clean_name, '-', '_');
clean_name = strrep(clean_name, '%', 'pct');
clean_name = strrep(clean_name, '(', '_');
clean_name = strrep(clean_name, ')', '_');
clean_name = strrep(clean_name, ' ', '_');
clean_name = strrep(clean_name, '.', '_');
clean_name = regexprep(clean_name, '[^a-zA-Z0-9_]', '_');

% Ensure it starts with a letter
if ~isempty(clean_name) && ~isletter(clean_name(1))
    clean_name = ['mat_' clean_name];
end

% Ensure it is not empty
if isempty(clean_name)
    clean_name = 'mat_unknown';
end
end


function wait_for_parallel_pool()
% WAIT_FOR_PARALLEL_POOL - Check parallel pool status (simplified version)
% Since pool is initialized synchronously in GUI, this function mainly reports status

% Check if pool exists
pool = gcp('nocreate');
if ~isempty(pool)
    fprintf('✓ Parallel pool ready (%d workers)\n', pool.NumWorkers);
    return;
end

% Check if Parallel Computing Toolbox is available
if ~license('test', 'Distrib_Computing_Toolbox')
    fprintf('⚠ Parallel Computing Toolbox not available - using serial computation\n');
    return;
end

% No pool found - try to create one immediately as fallback
fprintf('⚠ No parallel pool found - attempting to create fallback pool\n');
try
    num_physical_cores = feature('numcores');
    optimal_workers = max(2, num_physical_cores - 1);
    pool = parpool('local', optimal_workers);
    fprintf('✓ Fallback parallel pool created (%d workers)\n', pool.NumWorkers);
catch ME
    fprintf('⚠ Failed to create fallback parallel pool: %s\n', ME.message);
    fprintf('  Using serial computation\n');
end
end







