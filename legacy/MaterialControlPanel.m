function MaterialControlPanel()
%% MATERIALCONTROLPANEL - Advanced GUI for LGK Dendritic Growth Calculations
% ========================================================================
%
% FUNCTION PURPOSE:
% This is the main graphical user interface (GUI) for the LGK dendritic 
% growth model. It provides a comprehensive control panel for material 
% parameter selection, calculation parameter configuration, and computational 
% management with real-time progress monitoring.
%
% GUI ARCHITECTURE:
% The interface is organised into six main panels:
% 1. Material Selection - Choose from predefined materials or define custom
% 2. Material Parameters - Display/edit thermophysical properties
% 3. Calculation Parameters - Configure undercooling ranges and sequences
% 4. Initial Guess Functions - Set parametric or analytical approximations
% 5. Control Panel - Start/stop calculations and save results
% 6. Progress Display - Real-time calculation progress with time estimates
%
% SUPPORTED MATERIALS:
% - AZ91 (Magnesium alloy)
% - Al-4wt%Cu, Al-Cu, Al-Fe (Aluminium alloys)
% - Sn-Ag, Sn-Cu (Tin-based solders)
% - SCN-Acetone (Transparent model system from LGK experiments)
% - Mg-Alloy (General magnesium systems)
% - Custom (User-defined parameters)
%
% CALCULATION FEATURES:
% - Parallel processing with automatic worker pool management
% - Adaptive deltaT sequence generation with density control
% - Multiple initial guess strategies (parametric vs analytical)
% - Real-time progress monitoring with time estimation
% - Comprehensive result export to Excel and figure formats
%
% PARAMETER PERSISTENCE:
% - Automatically saves and restores material-specific parameters
% - Maintains separate parameter sets for each material
% - Preserves calculation settings between sessions
%
% INPUTS: None (GUI-driven parameter input)
% OUTPUTS: None (results saved to base workspace and files)
%
% DEPENDENCIES:
% - calculate_VR.m: Core Newton-Raphson solver
% - run_main_calculation.m: Multi-concentration calculation orchestrator
% - create_deltaT_sequence.m: Undercooling sequence generator
% - create_analytical_guess_functions.m: Analytical approximation functions
% - cleanup_parallel_pool.m: Parallel computing resource management
%
% USAGE:
% Simply call MaterialControlPanel() to launch the interface
% All interactions are handled through the GUI components
% ========================================================================

% Create main figure with reduced height
fig = figure(999);
set(fig, 'Position', [300, 200, 900, 650], ...
    'Name', 'Material Parameter Control Panel', ...
    'NumberTitle', 'off', ...
    'Tag', 'MaterialControlPanel', ...
    'MenuBar', 'none', ...
    'Toolbar', 'none', ...
    'Color', [0.94, 0.94, 0.94], ...
    'CloseRequestFcn', @close_panel);

% Store handles in figure's UserData
handles = struct();
handles.figure = fig;

% Create GUI components
handles = create_gui_components(fig, handles);

% Store handles in figure first
set(fig, 'UserData', handles);

fprintf('=== MAIN FUNCTION DEBUG ===\n');
fprintf('About to call load_saved_parameters\n');

% Load saved parameters if they exist, otherwise use defaults
loaded_params = load_saved_parameters(handles);

if loaded_params
    fprintf('Previous session parameters loaded successfully\n');
else
    % No saved parameters found, use default material with analytical approximation enabled
    set(handles.material_popup, 'Value', 1);
    set(handles.use_analytical_checkbox, 'Value', 1);
    fprintf('Using default material parameters\n');
end

% Apply the material selection (this will load material-specific defaults)
material_popup_Callback(handles.material_popup, [], handles);

% Apply analytical approximation setting
analytical_checkbox_Callback(handles.use_analytical_checkbox, [], handles);

% Ensure split point validation is applied
split_point_edit_Callback(handles.split_pt_edit, [], handles);

end


function handles = create_gui_components(fig, handles)
%% CREATE_GUI_COMPONENTS - Construct all GUI interface elements
%
% PURPOSE:
% Creates and positions all GUI panels, controls, and display elements
% with unified height variables for consistent layout management.
%
% PANEL STRUCTURE:
% - Material Selection (blue): Material type and thermal undercooling options
% - Material Parameters (yellow): Thermophysical property display/editing
% - Calculation Parameters (green): Undercooling ranges and sequence control
% - Initial Guess Functions (purple): Parametric function coefficients
% - Control Panel (grey): Calculation execution and file management
% - Progress Display (light grey): Real-time calculation progress
%
% INPUTS:
%   fig - Parent figure handle
%   handles - Structure for storing GUI component handles
%
% OUTPUTS:
%   handles - Updated structure containing all GUI component handles

% Material Selection Panel heights
mat_select_y = 10;  % Material selection row height

% Material Parameters Panel heights
param_row1_y = 110;  % First parameter row
param_row2_y = 80;   % Second parameter row
param_row3_y = 50;   % Third parameter row
param_row4_y = 20;   % Fourth parameter row

% Calculation Parameters Panel heights
calc_row1_y = 70;   % First calculation row
calc_row2_y = 40;    % Second calculation row
calc_c0_y = 10;      % C0 values row

% Initial Guess Functions Panel heights
v_function_y = 45;   % V function row
r_function_y = 15;   % R function row
auto_btn_y = 35;     % Auto button row

% Control Panel height
control_btn_y = 20;     % Control buttons row
instruction_y = 0;      % Instruction text row

% Progress Panel heights
progress_display_y = 15; % Progress display row

% Material Selection Panel
material_panel = uipanel('Parent', fig, ...
    'Title', 'Material Selection', ...
    'FontSize', 12, 'FontWeight', 'bold', ...
    'Position', [0.02, 0.89, 0.96, 0.09], ...
    'BackgroundColor', [0.9, 0.95, 1.0], ...
    'BorderType', 'line', 'BorderWidth', 2);

uicontrol('Parent', material_panel, ...
    'Style', 'text', ...
    'String', 'Select Material:', ...
    'Position', [10, mat_select_y-2, 120, 25], ...
    'FontSize', 11, 'FontWeight', 'bold', ...
    'BackgroundColor', [0.9, 0.95, 1.0]);

handles.material_popup = uicontrol('Parent', material_panel, ...
    'Style', 'popup', ...
    'String', {'AZ91', 'Al-4wt%Cu', 'Al-Cu', 'Al-Fe', 'Sn-Ag', 'Sn-Cu', 'SCN-Acetone(LGK 1984)', 'Mg-Alloy(Lin 2009)', 'Custom'}, ...
    'Position', [140, mat_select_y, 150, 25], ...
    'FontSize', 10, ...
    'Callback', {@material_popup_Callback});

% Add thermal undercooling option
uicontrol('Parent', material_panel, ...
    'Style', 'text', ...
    'String', 'Include Thermal Undercooling', ...
    'Position', [425, mat_select_y-3, 300, 25], ...
    'FontSize', 10, 'FontWeight', 'bold', ...
    'BackgroundColor', [0.9, 0.95, 1.0], ...
    'HorizontalAlignment', 'left');

handles.thermal_undercooling_checkbox = uicontrol('Parent', material_panel, ...
    'Style', 'checkbox', ...
    'String', '', ...
    'Position', [400, mat_select_y, 25, 25], ...
    'FontSize', 10, ...
    'Value', 1, ...
    'BackgroundColor', [0.9, 0.95, 1.0], ...
    'Callback', {@thermal_undercooling_Callback});

% Material Parameters Panel
param_panel = uipanel('Parent', fig, ...
    'Title', 'Material Parameters', ...
    'FontSize', 12, 'FontWeight', 'bold', ...
    'Position', [0.02, 0.64, 0.96, 0.24], ...
    'BackgroundColor', [1.0, 0.95, 0.9], ...
    'BorderType', 'line', 'BorderWidth', 2);

% Create parameter displays with full names, symbols, and units
param_info = {
    struct('name', 'Specific Heat Capacity', 'symbol', 'Cp', 'unit_vol', '(J/m³·K)', 'unit_spec', '(J/kg·K)', 'tag', 'cp_edit'), ...
    struct('name', 'Latent Heat', 'symbol', 'ΔH', 'unit_vol', '(J/m³)', 'unit_spec', '(J/kg)', 'tag', 'dh_edit'), ...
    struct('name', 'Thermal Diffusivity', 'symbol', 'a', 'unit_vol', '(m²/s)', 'unit_spec', '(m²/s)', 'tag', 'a_edit'), ...
    struct('name', 'Gibbs-Thomson Coefficient', 'symbol', 'Γ', 'unit_vol', '(K·m)', 'unit_spec', '(K·m)', 'tag', 'gibbs_edit'), ...
    struct('name', 'Diffusion Coefficient', 'symbol', 'D', 'unit_vol', '(m²/s)', 'unit_spec', '(m²/s)', 'tag', 'dal_edit'), ...
    struct('name', 'Liquidus Slope', 'symbol', 'm', 'unit_vol', '(K/wt%)', 'unit_spec', '(K/wt%)', 'tag', 'mal_edit'), ...
    struct('name', 'Partition Coefficient', 'symbol', 'k₀', 'unit_vol', '(-)', 'unit_spec', '(-)', 'tag', 'k0al_edit'), ...
    struct('name', 'Stability Constant', 'symbol', 'σ', 'unit_vol', '(-)', 'unit_spec', '(-)', 'tag', 'sigma_edit')
};

% Store parameter info for later use
handles.param_info = param_info;

% Create parameter controls in a 2-column layout using unified height variables
param_row_heights = [param_row1_y, param_row2_y, param_row3_y, param_row4_y];
for i = 1:length(param_info)
    col = mod(i-1, 2) + 1;  % 1 or 2
    row = ceil(i/2);        % 1, 2, 3, 4
    
    % Use unified height variable
    x_start = 20 + (col-1) * 420;  % Left column at 20, right at 440
    y_start = param_row_heights(row);
    
    % Create full name label (first line)
    handles.([param_info{i}.tag '_name']) = uicontrol('Parent', param_panel, ...
        'Style', 'text', ...
        'String', param_info{i}.name, ...
        'Position', [x_start, y_start, 200, 18], ...
        'FontSize', 9, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', [1.0, 0.95, 0.9]);
    
    % Create symbol and unit label (second line)
    unit_str = param_info{i}.unit_vol;  % Default to volumetric units
    symbol_unit_str = [param_info{i}.symbol, ' ', unit_str];
    
    handles.([param_info{i}.tag '_symbol']) = uicontrol('Parent', param_panel, ...
        'Style', 'text', ...
        'String', symbol_unit_str, ...
        'Position', [x_start, y_start-12, 200, 16], ...
        'FontSize', 8, ...
        'HorizontalAlignment', 'left', ...
        'ForegroundColor', [0.4, 0.4, 0.7], ...
        'BackgroundColor', [1.0, 0.95, 0.9]);
    
    % Create edit box (centered vertically between the two text lines)
    handles.(param_info{i}.tag) = uicontrol('Parent', param_panel, ...
        'Style', 'edit', ...
        'Position', [x_start + 210, y_start-6, 120, 22], ...
        'FontSize', 9, ...
        'BackgroundColor', 'white');
end

calc_panel = uipanel('Parent', fig, ...
    'Title', 'Calculation Parameters', ...
    'FontSize', 12, 'FontWeight', 'bold', ...
    'Position', [0.02, 0.44, 0.96, 0.19], ... 
    'BackgroundColor', [0.95, 1.0, 0.9], ...
    'BorderType', 'line', 'BorderWidth', 2);

% First row: ΔT_min, ΔT_max, sampling_interval
row1_names = {'ΔT_min:', 'ΔT_max:', 'Sampling Interval:'};
row1_tags = {'dt_min_edit', 'dt_max_edit', 'sampling_interval_edit'};
row1_defaults = {'0.1', '20', '0.5'};

for i = 1:length(row1_names)
    x_pos = 20 + (i-1) * 280;
    y_pos = calc_row1_y;
    
    uicontrol('Parent', calc_panel, ...
        'Style', 'text', ...
        'String', row1_names{i}, ...
        'Position', [x_pos, y_pos, 130, 20], ...
        'FontSize', 9, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', [0.95, 1.0, 0.9]);
    
    handles.(row1_tags{i}) = uicontrol('Parent', calc_panel, ...
        'Style', 'edit', ...
        'String', row1_defaults{i}, ...
        'Position', [x_pos + 135, y_pos, 100, 22], ...
        'FontSize', 9, ...
        'BackgroundColor', 'white');
    
    % Add callbacks for deltaT range validation
    if strcmp(row1_tags{i}, 'dt_min_edit') || strcmp(row1_tags{i}, 'dt_max_edit')
        set(handles.(row1_tags{i}), 'Callback', {@deltaT_range_edit_Callback});
    end
end

% Second row: split_point, lower_density_mult, upper_density_mult
row2_names = {'Split Point:', 'Lower Density Mult:', 'Upper Density Mult:'};
row2_tags = {'split_pt_edit', 'lower_density_mult_edit', 'upper_density_mult_edit'};
row2_defaults = {'1.0', '5.0', '1.0'};

for i = 1:length(row2_names)
    x_pos = 20 + (i-1) * 280;
    y_pos = calc_row2_y;
    
    uicontrol('Parent', calc_panel, ...
        'Style', 'text', ...
        'String', row2_names{i}, ...
        'Position', [x_pos, y_pos, 130, 20], ...
        'FontSize', 9, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', [0.95, 1.0, 0.9]);
    
    handles.(row2_tags{i}) = uicontrol('Parent', calc_panel, ...
        'Style', 'edit', ...
        'String', row2_defaults{i}, ...
        'Position', [x_pos + 135, y_pos, 100, 22], ...
        'FontSize', 9, ...
        'BackgroundColor', 'white');
    
    % Add callbacks for validation
    if strcmp(row2_tags{i}, 'split_pt_edit')
        set(handles.(row2_tags{i}), 'Callback', {@split_point_edit_Callback});
    end
end

% Third row: low_end_density, high_end_density, C0_values
row3_names = {'Low End Density:', 'High End Density:', 'C0 values:'};
row3_tags = {'low_density_edit', 'high_density_edit', 'c0_values_edit'};
row3_defaults = {'3.0', '1.0', '5:5:10'};

for i = 1:length(row3_names)
    x_pos = 20 + (i-1) * 280;
    y_pos = calc_c0_y;
    
    uicontrol('Parent', calc_panel, ...
        'Style', 'text', ...
        'String', row3_names{i}, ...
        'Position', [x_pos, y_pos, 130, 20], ...
        'FontSize', 9, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', [0.95, 1.0, 0.9]);
    
    handles.(row3_tags{i}) = uicontrol('Parent', calc_panel, ...
        'Style', 'edit', ...
        'String', row3_defaults{i}, ...
        'Position', [x_pos + 135, y_pos, 100, 22], ...
        'FontSize', 9, ...
        'BackgroundColor', 'white');
end

% Initial Guess Functions Panel
guess_panel = uipanel('Parent', fig, ...
    'Title', 'Initial Guess Functions', ...
    'FontSize', 12, 'FontWeight', 'bold', ...
    'Position', [0.02, 0.28, 0.96, 0.15], ...
    'BackgroundColor', [0.95, 0.9, 1.0], ...
    'BorderType', 'line', 'BorderWidth', 2);

% V guess function using unified height variable
uicontrol('Parent', guess_panel, ...
    'Style', 'text', ...
    'String', 'V_guess = (ΔT/(a×C0^b))^c', ...
    'Position', [10, v_function_y, 200, 20], ...
    'FontSize', 10, 'FontWeight', 'bold', ...
    'BackgroundColor', [0.95, 0.9, 1.0]);

uicontrol('Parent', guess_panel, ...
    'Style', 'text', ...
    'String', 'a:', ...
    'Position', [230, v_function_y, 20, 20], ...
    'FontSize', 9, ...
    'BackgroundColor', [0.95, 0.9, 1.0]);

handles.v_a_edit = uicontrol('Parent', guess_panel, ...
    'Style', 'edit', ...
    'String', '36.3', ...
    'Position', [255, v_function_y, 80, 22], ...
    'FontSize', 9, ...
    'BackgroundColor', 'white');

uicontrol('Parent', guess_panel, ...
    'Style', 'text', ...
    'String', 'b:', ...
    'Position', [345, v_function_y, 20, 20], ...
    'FontSize', 9, ...
    'BackgroundColor', [0.95, 0.9, 1.0]);

handles.v_b_edit = uicontrol('Parent', guess_panel, ...
    'Style', 'edit', ...
    'String', '0.695', ...
    'Position', [370, v_function_y, 80, 22], ...
    'FontSize', 9, ...
    'BackgroundColor', 'white');

uicontrol('Parent', guess_panel, ...
    'Style', 'text', ...
    'String', 'c:', ...
    'Position', [460, v_function_y, 20, 20], ...
    'FontSize', 9, ...
    'BackgroundColor', [0.95, 0.9, 1.0]);

handles.v_c_edit = uicontrol('Parent', guess_panel, ...
    'Style', 'edit', ...
    'String', '2.71', ...
    'Position', [485, v_function_y, 80, 22], ...
    'FontSize', 9, ...
    'BackgroundColor', 'white');

% R guess function using unified height variable
uicontrol('Parent', guess_panel, ...
    'Style', 'text', ...
    'String', 'R_guess = d × C0^e × ΔT^f', ...
    'Position', [10, r_function_y, 200, 20], ...
    'FontSize', 10, 'FontWeight', 'bold', ...
    'BackgroundColor', [0.95, 0.9, 1.0]);

uicontrol('Parent', guess_panel, ...
    'Style', 'text', ...
    'String', 'd:', ...
    'Position', [230, r_function_y, 20, 20], ...
    'FontSize', 9, ...
    'BackgroundColor', [0.95, 0.9, 1.0]);

handles.rd_edit = uicontrol('Parent', guess_panel, ...
    'Style', 'edit', ...
    'String', '3.93966e-5', ...
    'Position', [255, r_function_y, 80, 22], ...
    'FontSize', 9, ...
    'BackgroundColor', 'white');

uicontrol('Parent', guess_panel, ...
    'Style', 'text', ...
    'String', 'e:', ...
    'Position', [345, r_function_y, 20, 20], ...
    'FontSize', 9, ...
    'BackgroundColor', [0.95, 0.9, 1.0]);

handles.re_edit = uicontrol('Parent', guess_panel, ...
    'Style', 'edit', ...
    'String', '0', ...
    'Position', [370, r_function_y, 80, 22], ...
    'FontSize', 9, ...
    'BackgroundColor', 'white');

uicontrol('Parent', guess_panel, ...
    'Style', 'text', ...
    'String', 'f:', ...
    'Position', [460, r_function_y, 20, 20], ...
    'FontSize', 9, ...
    'BackgroundColor', [0.95, 0.9, 1.0]);

handles.rf_edit = uicontrol('Parent', guess_panel, ...
    'Style', 'edit', ...
    'String', '-1.2444', ...
    'Position', [485, r_function_y, 80, 22], ...
    'FontSize', 9, ...
    'BackgroundColor', 'white');

% Checkbox for using analytical approximation
handles.use_analytical_checkbox = uicontrol('Parent', guess_panel, ...
    'Style', 'checkbox', ...
    'String', 'Use Analytical Approximation', ...
    'Position', [630, auto_btn_y+5, 200, 20], ...
    'FontSize', 9, 'FontWeight', 'bold', ...
    'BackgroundColor', [0.95, 0.9, 1.0], ...
    'Value', 1, ...
    'Callback', {@analytical_checkbox_Callback});

% Help text for analytical approximation
uicontrol('Parent', guess_panel, ...
    'Style', 'text', ...
    'String', 'Uses Eqs. (8.91) & (8.92) from Solidification Book', ...
    'Position', [630, auto_btn_y-15, 200, 15], ...
    'FontSize', 7, ...
    'ForegroundColor', [0.5, 0.5, 0.5], ...
    'BackgroundColor', [0.95, 0.9, 1.0], ...
    'HorizontalAlignment', 'center');

% Control Buttons Panel
button_panel = uipanel('Parent', fig, ...
    'Title', 'Control Panel', ...
    'FontSize', 12, 'FontWeight', 'bold', ...
    'Position', [0.02, 0.14, 0.96, 0.13], ...
    'BackgroundColor', [0.9, 0.9, 0.9], ...
    'BorderType', 'line', 'BorderWidth', 2);

% Calculate button positions for centering (3 buttons now)
button_width = 130;
button_spacing = 40;
total_button_width = 3 * button_width + 2 * button_spacing + 20;
start_x = (900 - total_button_width) / 3;

% Control buttons using unified height variable
handles.start_btn = uicontrol('Parent', button_panel, ...
    'Style', 'pushbutton', ...
    'String', 'Start Calculation', ...
    'Position', [start_x, control_btn_y, button_width, 40], ...
    'FontSize', 10, 'FontWeight', 'bold', ...
    'BackgroundColor', [0.7, 1.0, 0.7], ...
    'ForegroundColor', [0, 0, 0], ...
    'Callback', {@start_calculation_Callback});

handles.stop_btn = uicontrol('Parent', button_panel, ...
    'Style', 'pushbutton', ...
    'String', 'Stop Calculation', ...
    'Position', [start_x + button_width + 10, control_btn_y, button_width, 40], ...
    'FontSize', 10, 'FontWeight', 'bold', ...
    'BackgroundColor', [1.0, 0.7, 0.7], ...
    'ForegroundColor', [0, 0, 0], ...
    'Enable', 'off', ...
    'Callback', {@stop_calculation_Callback});

handles.save_btn = uicontrol('Parent', button_panel, ...
    'Style', 'pushbutton', ...
    'String', 'Save Figures & Data', ...
    'Position', [start_x + 2*(button_width + 10), control_btn_y, button_width+20, 40], ...
    'FontSize', 10, 'FontWeight', 'bold', ...
    'BackgroundColor', [0.7, 0.9, 1.0], ...
    'ForegroundColor', [0, 0, 0], ...
    'Callback', {@save_data_Callback});

handles.close_btn = uicontrol('Parent', button_panel, ...
    'Style', 'pushbutton', ...
    'String', 'Close All Figures', ...
    'Position', [start_x + 3*(button_width + 10) + 20, control_btn_y, button_width, 40], ...
    'FontSize', 10, 'FontWeight', 'bold', ...
    'BackgroundColor', [0.7, 0.7, 1.0], ...
    'ForegroundColor', [0, 0, 0], ...
    'Callback', {@close_all_Callback});

% Add instruction text at bottom using unified height variable
handles.instruction_text = uicontrol('Parent', button_panel, ...
    'Style', 'text', ...
    'String', '', ...
    'Position', [20, instruction_y, 860, 15], ...
    'FontSize', 8, ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', [0.9, 0.9, 0.9], ...
    'ForegroundColor', [0.5, 0.5, 0.5], ...
    'Visible', 'off');

% Progress Panel
progress_panel = uipanel('Parent', fig, ...
    'Title', 'Progress Display', ...
    'FontSize', 12, 'FontWeight', 'bold', ...
    'Position', [0.02, 0.02, 0.96, 0.11], ... 
    'BackgroundColor', [0.95, 0.95, 0.95], ...
    'BorderType', 'line', 'BorderWidth', 2);

uicontrol('Parent', progress_panel, ...
    'Style', 'text', ...
    'String', 'Progress:', ...
    'Position', [20, progress_display_y, 80, 20], ...
    'FontSize', 10, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', [0.95, 0.95, 0.95]);

handles.progress_text_bar = uicontrol('Parent', progress_panel, ...
    'Style', 'text', ...
    'String', '[░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░]', ...
    'Position', [110, progress_display_y+2, 420, 16], ...
    'FontSize', 8, 'FontWeight', 'bold', ...
    'FontName', 'Courier New', ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', [0.95, 0.95, 0.95], ...
    'ForegroundColor', [0.2, 0.8, 0.2]);

handles.progress_percent = uicontrol('Parent', progress_panel, ...
    'Style', 'text', ...
    'String', '0%', ...
    'Position', [535, progress_display_y, 60, 20], ...
    'FontSize', 9, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', [0.95, 0.95, 0.95]);

handles.progress_elapsed = uicontrol('Parent', progress_panel, ...
    'Style', 'text', ...
    'String', 'Elapsed: 0.00 min', ...
    'Position', [605, progress_display_y, 120, 20], ...
    'FontSize', 9, ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', [0.95, 0.95, 0.95]);

handles.progress_remaining = uicontrol('Parent', progress_panel, ...
    'Style', 'text', ...
    'String', 'Remaining: -- min', ...
    'Position', [725, progress_display_y, 150, 20], ...
    'FontSize', 9, ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', [0.95, 0.95, 0.95]);

% Add author information in top-right corner (very subtle)
uicontrol('Parent', fig, ...
    'Style', 'text', ...
    'String', 'Author: Xinyi Hao', ...
    'Position', [805, 620, 75, 15], ...
    'FontSize', 7, ...
    'FontWeight', 'normal', ...
    'HorizontalAlignment', 'right', ...
    'BackgroundColor', [0.94, 0.94, 0.94], ...
    'ForegroundColor', [0.6, 0.6, 0.6], ...
    'Enable', 'inactive');

% Initial validation of split point
split_point_edit_Callback(handles.split_pt_edit, [], handles);

% Initial validation of split point
split_point_edit_Callback(handles.split_pt_edit, [], handles);

end






function material_popup_Callback(hObject, eventdata, ~)
%% MATERIAL_POPUP_CALLBACK - Handle material selection changes
%
% PURPOSE:
% Updates all material parameters, enables/disables appropriate controls,
% and loads material-specific calculation settings when user changes
% the selected material type.
%
% FUNCTIONALITY:
% - Loads predefined thermophysical parameters for selected material
% - Updates parameter labels with appropriate units (volumetric/specific)
% - Enables/disables parameter editing based on material type
% - Applies material-specific initial guess function coefficients
% - Handles special cases (SCN-Acetone polynomial interpolation)
%
% MATERIAL-SPECIFIC BEHAVIOUR:
% - Predefined materials: Parameters locked, calculation settings enabled
% - SCN-Acetone: Uses polynomial interpolation, limited deltaT range
% - Custom material: All parameters and settings editable

fig = ancestor(hObject, 'figure');
handles = get(fig, 'UserData');

% Load saved parameters for the newly selected material
if exist('material_memory.mat', 'file')
    load('material_memory.mat', 'material_params');
    selection = get(hObject, 'Value');
    load_material_specific_parameters(handles, selection, material_params);
end

materials = get_material_parameters();
selection = get(hObject, 'Value');
material_names = {'AZ91', 'Al_4wtCu', 'Al_Cu', 'Al_Fe', 'Sn_Ag', 'Sn_Cu', 'SCN_Acetone', 'Mg_Alloy', 'Custom'};

if selection <= 8
    material = materials.(material_names{selection});
    
    % Update parameter fields (removed DSv)
    set(handles.cp_edit, 'String', num2str(material.Cp, '%.6e'));
    set(handles.dh_edit, 'String', num2str(material.DH, '%.6e'));
    set(handles.a_edit, 'String', num2str(material.a, '%.6e'));
    set(handles.gibbs_edit, 'String', num2str(material.Gibbs_Tom, '%.6e'));
    set(handles.dal_edit, 'String', num2str(material.D_Al, '%.6e'));
    set(handles.mal_edit, 'String', num2str(material.m_Al, '%.6f'));
    set(handles.k0al_edit, 'String', num2str(material.k0_Al, '%.6f'));
    set(handles.sigma_edit, 'String', num2str(material.sigma, '%.6f'));
        
    % Update V_guess parameters based on material type
    if selection == 1  % AZ91
        set(handles.v_a_edit, 'String', '43.7305');
        set(handles.v_b_edit, 'String', '0.640395');
        set(handles.v_c_edit, 'String', '2.6745');
        % Set R_guess parameters for Al-based alloys
        set(handles.rd_edit, 'String', '1.74591e-5');
        set(handles.re_edit, 'String', '0.360984');
        set(handles.rf_edit, 'String', '-1.36146');
    elseif selection == 2 || selection == 3 || selection == 4 %  Al-4wt%Cu, Al-Cu, Al-Fe
        set(handles.v_a_edit, 'String', '52.8963785');% =(1/49.14)^(1/2.5)*(10^6)^(1/2.5)
        set(handles.v_b_edit, 'String', '0.6');
        set(handles.v_c_edit, 'String', '2.5');
        % Set R_guess parameters for Sn-based alloys
        set(handles.rd_edit, 'String', '20.34e-6');
        set(handles.re_edit, 'String', '0.25');
        set(handles.rf_edit, 'String', '-1.25');
    elseif selection == 5 || selection == 6  % Sn-Ag, Sn-Cu (tin-based)
        set(handles.v_a_edit, 'String', '27.0497');
        set(handles.v_b_edit, 'String', '0.630488');
        set(handles.v_c_edit, 'String', '2.85702');
        % Set R_guess parameters for Sn-based alloys
        set(handles.rd_edit, 'String', '5.21454e-06');
        set(handles.re_edit, 'String', '0.62773');
        set(handles.rf_edit, 'String', '-1.60336');
    elseif selection == 7  % SCN-Acetone
        % For SCN-Acetone, disable the V_guess and R_guess parameter editing
        % as it uses fitted polynomial and power function interpolation
        set(handles.v_a_edit, 'String', 'Fitted');
        set(handles.v_b_edit, 'String', 'Polynomial');
        set(handles.v_c_edit, 'String', 'Degree 6');
        set(handles.rd_edit, 'String', 'Fitted');
        set(handles.re_edit, 'String', 'Power');
        set(handles.rf_edit, 'String', 'Function');
        
        % Disable editing for polynomial mode
        set(handles.v_a_edit, 'Enable', 'off');
        set(handles.v_b_edit, 'Enable', 'off');
        set(handles.v_c_edit, 'Enable', 'off');
        set(handles.rd_edit, 'Enable', 'off');
        set(handles.re_edit, 'Enable', 'off');
        set(handles.rf_edit, 'Enable', 'off');
        
        % Set fixed calculation parameters for SCN-Acetone
        set(handles.dt_min_edit, 'String', '0.5');
        set(handles.dt_max_edit, 'String', '0.9');
        set(handles.sampling_interval_edit, 'String', '0.4');
        
        % Disable calculation parameter editing for SCN-Acetone (except C0 values)
        set(handles.dt_min_edit, 'Enable', 'off');
        set(handles.dt_max_edit, 'Enable', 'off');
        set(handles.sampling_interval_edit, 'Enable', 'off');
        set(handles.split_pt_edit, 'Enable', 'off');
        set(handles.lower_density_mult_edit, 'Enable', 'off');
        set(handles.upper_density_mult_edit, 'Enable', 'off');
        set(handles.low_density_edit, 'Enable', 'off');
        set(handles.high_density_edit, 'Enable', 'off');
        % Keep C0 values enabled
        set(handles.c0_values_edit, 'Enable', 'on');
    elseif selection == 8  % Mg-Alloy
        set(handles.v_a_edit, 'String', '34.436866');    
        set(handles.v_b_edit, 'String', '0.7066862');    
        set(handles.v_c_edit, 'String', '2.7000614');    
        % Set R_guess parameters for Mg-alloy
        set(handles.rd_edit, 'String', '1.0985739e-05');
        set(handles.re_edit, 'String', '0.52125156');
        set(handles.rf_edit, 'String', '-1.3846376');
    end
    
    % Update the labels with appropriate units
    update_parameter_labels(handles, material.unit_type);
    
    % Disable editing for predefined materials (removed DSv)
    fields = {'cp_edit', 'dh_edit', 'a_edit', 'gibbs_edit', ...
              'dal_edit', 'mal_edit', 'k0al_edit', 'sigma_edit'};
    for i = 1:length(fields)
        set(handles.(fields{i}), 'Enable', 'off');
    end
    
    % Keep V_guess and R_guess parameters always enabled for user modification (except for SCN-Acetone)
    if selection ~= 7  % Not SCN-Acetone
        % Enable V_guess parameters for non-SCN materials
        v_guess_fields = {'v_a_edit', 'v_b_edit', 'v_c_edit'};
        for i = 1:length(v_guess_fields)
            set(handles.(v_guess_fields{i}), 'Enable', 'on');
        end
        
        % Enable R_guess parameters for non-SCN materials
        r_guess_fields = {'rd_edit', 're_edit', 'rf_edit'};
        for i = 1:length(r_guess_fields)
            set(handles.(r_guess_fields{i}), 'Enable', 'on');
        end
    else
        % For SCN-Acetone, keep both V_guess and R_guess parameters disabled as they use polynomial interpolation
        v_guess_fields = {'v_a_edit', 'v_b_edit', 'v_c_edit'};
        for i = 1:length(v_guess_fields)
            set(handles.(v_guess_fields{i}), 'Enable', 'off');
        end
        
        r_guess_fields = {'rd_edit', 're_edit', 'rf_edit'};
        for i = 1:length(r_guess_fields)
            set(handles.(r_guess_fields{i}), 'Enable', 'off');
        end
    end
elseif selection == 9  % Custom material
    % Enable editing for custom material (removed DSv)
    fields = {'cp_edit', 'dh_edit', 'a_edit', 'gibbs_edit', ...
              'dal_edit', 'mal_edit', 'k0al_edit', 'sigma_edit'};
    for i = 1:length(fields)
        set(handles.(fields{i}), 'Enable', 'on');
    end
    
    % Enable V_guess parameter editing for custom material
    v_guess_fields = {'v_a_edit', 'v_b_edit', 'v_c_edit'};
    for i = 1:length(v_guess_fields)
        set(handles.(v_guess_fields{i}), 'Enable', 'on');
    end
    
    % Enable R_guess parameter editing for custom material
    r_guess_fields = {'rd_edit', 're_edit', 'rf_edit'};
    for i = 1:length(r_guess_fields)
        set(handles.(r_guess_fields{i}), 'Enable', 'on');
    end
    
    % Enable all calculation parameters for custom material
    set(handles.dt_min_edit, 'Enable', 'on');
    set(handles.dt_max_edit, 'Enable', 'on');
    set(handles.sampling_interval_edit, 'Enable', 'on');
    set(handles.split_pt_edit, 'Enable', 'on');
    set(handles.lower_density_mult_edit, 'Enable', 'on');
    set(handles.upper_density_mult_edit, 'Enable', 'on');
    set(handles.low_density_edit, 'Enable', 'on');
    set(handles.high_density_edit, 'Enable', 'on');
    set(handles.c0_values_edit, 'Enable', 'on');
    
    % Set default unit labels for custom material
    update_parameter_labels(handles, 'volumetric');
end

% For non-SCN materials, ensure all calculation parameters are enabled
if selection ~= 7 && selection ~= 9 % Not SCN and not Custom
    set(handles.dt_min_edit, 'Enable', 'on');
    set(handles.dt_max_edit, 'Enable', 'on');
    set(handles.sampling_interval_edit, 'Enable', 'on');
    set(handles.split_pt_edit, 'Enable', 'on');
    set(handles.lower_density_mult_edit, 'Enable', 'on');
    set(handles.upper_density_mult_edit, 'Enable', 'on');
    set(handles.low_density_edit, 'Enable', 'on');
    set(handles.high_density_edit, 'Enable', 'on');
    set(handles.c0_values_edit, 'Enable', 'on');
end

% Update figure's UserData
set(fig, 'UserData', handles);

end

function update_parameter_labels(handles, unit_type)
% Update parameter labels based on unit type

param_info = handles.param_info;

for i = 1:length(param_info)
    % Select appropriate unit based on type
    if strcmp(unit_type, 'volumetric')
        unit_str = param_info{i}.unit_vol;
    else
        unit_str = param_info{i}.unit_spec;
    end
    
    % Update the symbol and unit label
    symbol_unit_str = [param_info{i}.symbol, ' ', unit_str];
    set(handles.([param_info{i}.tag '_symbol']), 'String', symbol_unit_str);
end

end

function materials = get_material_parameters()
%% GET_MATERIAL_PARAMETERS - Retrieve predefined material property database
%
% PURPOSE:
% Provides comprehensive thermophysical property database for all
% supported materials, with appropriate unit systems and literature values.
%
% PARAMETER DEFINITIONS:
% - Cp/Cpv: Specific/volumetric heat capacity [J/kg·K or J/m³·K]
% - DH/DHv: Specific/volumetric latent heat [J/kg or J/m³]
% - a: Thermal diffusivity [m²/s]
% - Gibbs_Tom: Gibbs-Thomson coefficient [K·m]
% - D_Al: Interdiffusion coefficient [m²/s]
% - m_Al: Liquidus slope [K/wt%]
% - k0_Al: Partition coefficient [-]
% - sigma: Stability constant [-] = 1/(4π²)
% - unit_type: Indicates 'volumetric' or 'specific' unit system

% AZ91 - volumetric units
materials.AZ91.Cp = 1.8E+06;    % Specific heat (J/m³·K)
materials.AZ91.DH = 6.55E+08;   % Latent heat (J/m³)
materials.AZ91.a = 100/1.8E+06; % Diffusion coefficient (m²/s) = lambda/Cpv
% materials.AZ91.DSv = 595883;    % Entropy difference (J/m³·K)
materials.AZ91.Gibbs_Tom = 2*0.065/595883; % Gibbs-Thomson coefficient
materials.AZ91.D_Al = 3E-9;     % Diffusion coefficient (m²/s)
materials.AZ91.m_Al = -5.75;    % Liquidus slope (K/wt%)
materials.AZ91.k0_Al = 0.31;    % Partition coefficient
materials.AZ91.sigma = 1/(4*(pi^2)); % stability constant
materials.AZ91.unit_type = 'volumetric';

% Al-4wt%Cu - specific units (J/kg·K and J/kg)
materials.Al_4wtCu.Cp = 1070;           % Specific heat capacity (J/kg·K)
materials.Al_4wtCu.DH = 3.97e5;         % Latent heat (J/kg)
materials.Al_4wtCu.a = 4.313e-5;        % Thermal diffusivity (m²/s)
% materials.Al_4wtCu.DSv = 3.97e5/331.24; % Approximate entropy difference
materials.Al_4wtCu.Gibbs_Tom = 2.4e-7;  % Gibbs-Thomson coefficient (K·m)
materials.Al_4wtCu.D_Al = 3e-9;         % Interdiffusion coefficient (m²/s)
materials.Al_4wtCu.m_Al = -3.37;        % Liquidus slope (K/wt%)
materials.Al_4wtCu.k0_Al = 0.17;        % Partition coefficient
materials.Al_4wtCu.sigma = 1/(4*(pi^2)); % stability constant
materials.Al_4wtCu.unit_type = 'specific';

% Al-Cu - volumetric units
materials.Al_Cu.Cp = 2.67e6;     % Specific heat capacity (J/m³·K)
materials.Al_Cu.DH = 971e6;      % Latent heat (J/m³)
materials.Al_Cu.a = 0.34e-4;     % Thermal diffusivity (m²/s)
% materials.Al_Cu.DSv = 971e6/600; % Approximate entropy difference (using approx melting point)
materials.Al_Cu.Gibbs_Tom = 2.4e-7; % Gibbs-Thomson coefficient (K·m)
materials.Al_Cu.D_Al = 3e-9;     % Interdiffusion coefficient (m²/s)
materials.Al_Cu.m_Al = -3.37;    % Liquidus slope (K/wt%)
materials.Al_Cu.k0_Al = 0.17;    % Partition coefficient
materials.Al_Cu.sigma = 1/(4*(pi^2)); % stability constant
materials.Al_Cu.unit_type = 'volumetric';

% Al-Fe - volumetric units
materials.Al_Fe.Cp = 2.67e6;     % Specific heat capacity (J/m³·K)
materials.Al_Fe.DH = 971e6;      % Latent heat (J/m³)
materials.Al_Fe.a = 0.34e-4;     % Thermal diffusivity (m²/s)
% materials.Al_Fe.DSv = 971e6/600; % Approximate entropy difference
materials.Al_Fe.Gibbs_Tom = 1e-7; % Gibbs-Thomson coefficient (K·m)
materials.Al_Fe.D_Al = 2e-9;     % Interdiffusion coefficient (m²/s)
materials.Al_Fe.m_Al = -3.7;     % Liquidus slope (K/wt%)
materials.Al_Fe.k0_Al = 0.038;   % Partition coefficient
materials.Al_Fe.sigma = 1/(4*(pi^2)); % stability constant
materials.Al_Fe.unit_type = 'volumetric';

% Sn-Ag - specific units (J/kg·K and J/kg)
materials.Sn_Ag.Cp = 249;        % Specific heat capacity (J/kg·K)
materials.Sn_Ag.DH = 61810.62;   % Latent heat (J/kg)
materials.Sn_Ag.a = 1.5e-5;      % Thermal diffusivity (m²/s)
% materials.Sn_Ag.DSv = 61810.62/232; % Entropy difference (using Sn melting point)
materials.Sn_Ag.Gibbs_Tom = 8.54e-8; % Gibbs-Thomson coefficient (K·m)
materials.Sn_Ag.D_Al = 1.82e-9;  % Interdiffusion coefficient (m²/s)
materials.Sn_Ag.m_Al = -3.14;    % Liquidus slope (K/wt%)
materials.Sn_Ag.k0_Al = 0.0191;  % Partition coefficient
materials.Sn_Ag.sigma = 1/(4*(pi^2)); % stability constant
materials.Sn_Ag.unit_type = 'specific';

% Sn-Cu - specific units (J/kg·K and J/kg)
materials.Sn_Cu.Cp = 223;        % Specific heat capacity (J/kg·K)
materials.Sn_Cu.DH = 59212.54;   % Latent heat (J/kg)
materials.Sn_Cu.a = 2.54e-5;     % Thermal diffusivity (m²/s)
% materials.Sn_Cu.DSv = 59212.54/232; % Entropy difference
materials.Sn_Cu.Gibbs_Tom = 7.85e-8; % Gibbs-Thomson coefficient (K·m)
materials.Sn_Cu.D_Al = 0.86e-9;  % Interdiffusion coefficient (m²/s)
materials.Sn_Cu.m_Al = -5.79;    % Liquidus slope (K/wt%)
materials.Sn_Cu.k0_Al = 0.0051;  % Partition coefficient
materials.Sn_Cu.sigma = 1/(4*(pi^2)); % stability constant
materials.Sn_Cu.unit_type = 'specific';

% SCN-Acetone (LGK 1984) - specific units (J/kg·K and J/kg)
materials.SCN_Acetone.Cp = 1937.5;          % Specific heat capacity (J/kg·K)
materials.SCN_Acetone.DH = 46.26e3;         % Latent heat (J/kg)
materials.SCN_Acetone.a = 1.14e-7;          % Thermal diffusivity (m²/s)
% materials.SCN_Acetone.DSv = 46.26e3/331.24; % Entropy difference (using melting point)
materials.SCN_Acetone.Gibbs_Tom = 6.62e-8;  % Gibbs-Thomson coefficient (K·m)
materials.SCN_Acetone.D_Al = 1.27e-9;       % Diffusion coefficient (m²/s)
materials.SCN_Acetone.m_Al = -2.16;         % Liquidus slope (K/mol%)
materials.SCN_Acetone.k0_Al = 0.103;        % Partition coefficient
materials.SCN_Acetone.sigma = 1/(4*(pi^2)); % Stability constant
materials.SCN_Acetone.unit_type = 'specific';

% Mg alloy (Yin 2009) - specific units (J/kg·K and J/kg)
materials.Mg_Alloy.Cp = 1200;               % Specific heat capacity (J/kg·K)
materials.Mg_Alloy.DH = 3.7e5;              % Latent heat (J/kg)
materials.Mg_Alloy.a = 80/(1650*1350);         % Thermal diffusivity (m²/s) - estimated
% materials.Mg_Alloy.DSv = 3.7e5/705;         % Entropy difference (using eutectic temp)
materials.Mg_Alloy.Gibbs_Tom = 2.0e-7;      % Gibbs-Thomson coefficient (K·m)
materials.Mg_Alloy.D_Al = 5.0e-9;             % Diffusion coefficient (m²/s) - estimated
materials.Mg_Alloy.m_Al = -5.75;            % Liquidus slope (K/wt%) - estimated
materials.Mg_Alloy.k0_Al = 0.31;            % Partition coefficient - estimated
materials.Mg_Alloy.sigma = 1/(4*(pi^2));    % Stability constant
materials.Mg_Alloy.unit_type = 'specific';


end

function start_calculation_Callback(hObject, eventdata, ~)
%% START_CALCULATION_CALLBACK - Initiate LGK model calculations
%
% PURPOSE:
% Orchestrates the complete calculation workflow including parallel pool
% initialization, parameter validation, calculation execution, and result
% processing with comprehensive error handling.
%
% CALCULATION WORKFLOW:
% 1. Initialize parallel computing pool for performance optimization
% 2. Extract and validate all GUI parameters
% 3. Configure material parameters structure for solver
% 4. Execute multi-concentration calculations via run_main_calculation()
% 5. Process and store results for plotting and Excel export
% 6. Generate comprehensive visualizations and surface fitting
%
% PARALLEL PROCESSING:
% - Automatically detects available CPU cores
% - Creates optimal worker pool for calculation acceleration
% - Maintains pool between calculations for efficiency
% - Handles graceful fallback to serial computation if needed
%
% ERROR HANDLING:
% - Validates all input parameters before calculation
% - Provides user feedback for parameter errors
% - Maintains GUI responsiveness during long calculations
% - Cleans up resources on calculation failure

fig = ancestor(hObject, 'figure');

% Check if GUI is still available
if ~isvalid(fig) || ~ishandle(fig) || ~isgraphics(fig)
    fprintf('Error: GUI has been closed or destroyed\n');
    return;
end

handles = get(fig, 'UserData');

% Check handles
if ~isstruct(handles) || ~isfield(handles, 'start_btn')
    fprintf('Error: GUI handles are invalid\n');
    return;
end

% Check start_btn
if ~isvalid(handles.start_btn) || ~ishandle(handles.start_btn)
    fprintf('Error: Start button has been destroyed\n');
    return;
end

handles = get(fig, 'UserData');

try
    % Save current parameters before starting calculation
    save_current_parameters(handles);
    
    % Record parallel pool initialization time
    pool_init_start_time = now;
    pool_existed_before = ~isempty(gcp('nocreate'));
    
    % Initialize parallel pool when starting calculation
    initialize_parallel_pool_for_calculation();
    
    % Calculate actual initialization time
    pool_init_end_time = now;
    actual_pool_init_time = (pool_init_end_time - pool_init_start_time) * 24 * 60; % in minutes
    
    % Store the timing information for the calculation
    if pool_existed_before
        fprintf('Parallel pool was already available - no initialization time (%.3f min)\n', actual_pool_init_time);
        setappdata(fig, 'actual_pool_init_time', 0); % No real initialization time
        setappdata(fig, 'pool_was_pre_existing', true);
    else
        fprintf('Parallel pool created during calculation start - initialization time: %.2f min\n', actual_pool_init_time);
        setappdata(fig, 'actual_pool_init_time', actual_pool_init_time);
        setappdata(fig, 'pool_was_pre_existing', false);
    end
    
    % Automatically close all other figures before starting calculation
    fprintf('Closing all existing figures before starting new calculation...\n');
    close_all_Callback(handles.close_btn, [], handles);
    

    % Initialize global stop flag
    global STOP_CALCULATION;
    STOP_CALCULATION = false;
    
    clear functions;
    
    % Disable the start button during calculation and show instruction
    set(handles.start_btn, 'Enable', 'off');
    set(handles.start_btn, 'String', 'Calculating...');
    set(handles.start_btn, 'BackgroundColor', [0.8, 0.8, 0.8]);
    
    % Enable stop button during calculation
    set(handles.stop_btn, 'Enable', 'on');

    % Show instruction text
    set(handles.instruction_text, 'String', 'To force stop: Press Ctrl+C in Command Window (will lose current calculation progress)');
    set(handles.instruction_text, 'Visible', 'on');
    
    % Force GUI update
    drawnow;
    
    
    % Save current parameters
    save_current_parameters(handles);
    
    % Get material parameters (removed DSv)
    Cpv = str2double(get(handles.cp_edit, 'String'));
    DHv = str2double(get(handles.dh_edit, 'String'));
    a = str2double(get(handles.a_edit, 'String'));
    Gibbs_Tom = str2double(get(handles.gibbs_edit, 'String'));
    D_Al = str2double(get(handles.dal_edit, 'String'));
    m_Al = str2double(get(handles.mal_edit, 'String'));
    k0_Al = str2double(get(handles.k0al_edit, 'String'));
    sigma = str2double(get(handles.sigma_edit, 'String'));
    
    % Get calculation parameters
    deltaT_min = str2double(get(handles.dt_min_edit, 'String'));
    deltaT_max = str2double(get(handles.dt_max_edit, 'String'));
    sampling_interval = str2double(get(handles.sampling_interval_edit, 'String'));
    split_point = str2double(get(handles.split_pt_edit, 'String'));
    lower_density_mult = str2double(get(handles.lower_density_mult_edit, 'String'));
    upper_density_mult = str2double(get(handles.upper_density_mult_edit, 'String'));
    low_end_density = str2double(get(handles.low_density_edit, 'String'));
    high_end_density = str2double(get(handles.high_density_edit, 'String'));
    
    % Get C0 values
    c0_str = get(handles.c0_values_edit, 'String');
    C0_values = eval(c0_str);
    
    
    % Get current material name for export
    material_names = {'AZ91', 'Al-4wt%Cu', 'Al-Cu', 'Al-Fe', 'Sn-Ag', 'Sn-Cu', 'SCN-Acetone(LGK 1984)', 'Mg-Alloy(Lin 2009)', 'Custom'};
    material_selection = get(handles.material_popup, 'Value');
    current_material_name = material_names{material_selection};

    % Get thermal undercooling setting
    include_thermal_undercooling = get(handles.thermal_undercooling_checkbox, 'Value');
    
    % Get analytical approximation setting
    use_analytical_approximation = get(handles.use_analytical_checkbox, 'Value');
    
    % Get NEW V_guess parameters (handle analytical approximation case) - MOVED HERE
    if use_analytical_approximation
        % When using analytical approximation, GUI shows text instead of numbers
        % Use material-specific default values for plotting purposes
        material_selection = get(handles.material_popup, 'Value');
        
        if material_selection == 2 || material_selection == 3 || material_selection == 4 % Al-based materials
            V_a = 52.8964; V_b = 0.6; V_c = 2.5;
            R_d = 20.34e-6; R_e = 0.25; R_f = -1.25;
        elseif material_selection == 5 || material_selection == 6 % Sn-based materials
            V_a = 27.0497; V_b = 0.630488; V_c = 2.85702;
            R_d = 5.21454e-06; R_e = 0.62773; R_f = -1.60336;
        elseif material_selection == 1 || material_selection == 8 % Mg-based materials
            V_a = 43.7305; V_b = 0.640395; V_c = 2.6745;
            R_d = 1.74591e-5; R_e = 0.360984; R_f = -1.36146;
        else % Default to Al-Cu values
            V_a = 52.8964; V_b = 0.6; V_c = 2.5;
            R_d = 20.34e-6; R_e = 0.25; R_f = -1.25;
        end
        
        fprintf('Using material-specific default parameters for plotting:\n');
        fprintf('  V_a=%.3f, V_b=%.3f, V_c=%.3f\n', V_a, V_b, V_c);
        fprintf('  R_d=%.3e, R_e=%.3f, R_f=%.3f\n', R_d, R_e, R_f);
    else
        % For parametric approximation, get values from GUI
        V_a = str2double(get(handles.v_a_edit, 'String'));
        V_b = str2double(get(handles.v_b_edit, 'String'));
        V_c = str2double(get(handles.v_c_edit, 'String'));
        
        R_d = str2double(get(handles.rd_edit, 'String'));
        R_e = str2double(get(handles.re_edit, 'String'));
        R_f = str2double(get(handles.rf_edit, 'String'));
    end
    
    % Run the main calculation
    run_main_calculation(Cpv, DHv, a, Gibbs_Tom, D_Al, m_Al, k0_Al, sigma, ...
                    deltaT_min, deltaT_max, sampling_interval, split_point, ...
                    lower_density_mult, upper_density_mult, low_end_density, ...
                    high_end_density, C0_values, V_a, V_b, V_c, R_d, R_e, R_f, ...
                    include_thermal_undercooling, use_analytical_approximation);
    
    % Update material name and calculation parameters in export data
    if evalin('base', 'exist(''export_data'', ''var'')')
        export_data = evalin('base', 'export_data');
        export_data.material_name = current_material_name;
        
        % Store calculation parameters from GUI
        export_data.calc_params = struct();
        export_data.calc_params.deltaT_min = str2double(get(handles.dt_min_edit, 'String'));
        export_data.calc_params.deltaT_max = str2double(get(handles.dt_max_edit, 'String'));
        export_data.calc_params.sampling_interval = str2double(get(handles.sampling_interval_edit, 'String'));
        export_data.calc_params.split_point = str2double(get(handles.split_pt_edit, 'String'));
        export_data.calc_params.lower_density_mult = str2double(get(handles.lower_density_mult_edit, 'String'));
        export_data.calc_params.upper_density_mult = str2double(get(handles.upper_density_mult_edit, 'String'));
        export_data.calc_params.low_end_density = str2double(get(handles.low_density_edit, 'String'));
        export_data.calc_params.high_end_density = str2double(get(handles.high_density_edit, 'String'));
        export_data.calc_params.c0_values = get(handles.c0_values_edit, 'String');
        export_data.calc_params.V_a = str2double(get(handles.v_a_edit, 'String'));
        export_data.calc_params.V_b = str2double(get(handles.v_b_edit, 'String'));
        export_data.calc_params.V_c = str2double(get(handles.v_c_edit, 'String'));
        export_data.calc_params.R_d = str2double(get(handles.rd_edit, 'String'));
        export_data.calc_params.R_e = str2double(get(handles.re_edit, 'String'));
        export_data.calc_params.R_f = str2double(get(handles.rf_edit, 'String'));
        
        assignin('base', 'export_data', export_data);
    end
    
    % Re-enable the start button after calculation and hide instruction
    if isvalid(fig) && ishandle(fig) && isgraphics(fig)
        handles = get(fig, 'UserData');
        if isstruct(handles) && isfield(handles, 'start_btn') && ...
           isvalid(handles.start_btn) && ishandle(handles.start_btn)
            set(handles.start_btn, 'Enable', 'on');
            set(handles.start_btn, 'String', 'Start Calculation');
            set(handles.start_btn, 'BackgroundColor', [0.7, 1.0, 0.7]);
            
            % Disable stop button
            set(handles.stop_btn, 'Enable', 'off');

            % Hide instruction text
            set(handles.instruction_text, 'Visible', 'off');
        end
    end

    
    % Show success message
    msgbox('Calculation completed successfully! You can now save the results or run another calculation.', 'Success', 'modal');
    
catch ME
    % Print detailed error information
    fprintf('Error in calculation: %s\n', ME.message);
    fprintf('Error occurred in: %s at line %d\n', ME.stack(1).file, ME.stack(1).line);
    
    % Stop and cleanup progress timer if it exists
    try
        % Try to find and stop any running progress timer
        all_timers = timerfindall();
        for i = 1:length(all_timers)
            if isvalid(all_timers(i))
                stop(all_timers(i));
                delete(all_timers(i));
            end
        end
        fprintf('All timers cleaned up due to error.\n');
    catch
        % Silent cleanup failure
    end
    
    % Re-enable the start button if error occurs and hide instruction
    if isvalid(fig) && ishandle(fig) && isgraphics(fig)
        handles = get(fig, 'UserData');
        
        if isstruct(handles) && isfield(handles, 'start_btn') && ...
           isvalid(handles.start_btn) && ishandle(handles.start_btn)
            set(handles.start_btn, 'Enable', 'on');
            set(handles.start_btn, 'String', 'Start Calculation');
            set(handles.start_btn, 'BackgroundColor', [0.7, 1.0, 0.7]);
            
            % Hide instruction text
            set(handles.instruction_text, 'Visible', 'off');
        end
    end
    
    % Clean up parallel pool due to error
    try
        cleanup_parallel_pool();
        fprintf('Parallel pool cleaned up due to error.\n');
    catch
        fprintf('Could not clean up parallel pool.\n');
    end

    msgbox(['Error during calculation: ' ME.message], 'Error', 'error', 'modal');
end

end

function save_data_Callback(hObject, eventdata, ~)
% Save all figures and data including Excel export
% Get handles from figure
fig = ancestor(hObject, 'figure');
handles = get(fig, 'UserData');

try
    % Save current parameters before saving data
    save_current_parameters(handles);
    
    % Disable the save button and show the save status
    set(handles.save_btn, 'Enable', 'off');
    set(handles.save_btn, 'String', 'Saving...');
    set(handles.save_btn, 'BackgroundColor', [0.8, 0.8, 0.8]);
    
    % Force GUI update
    drawnow;
    
    % Select save directory
    save_dir = uigetdir(pwd, 'Select directory to save results');
    if save_dir == 0
        % If user cancels the process, restore the button state
        set(handles.save_btn, 'Enable', 'on');
        set(handles.save_btn, 'String', 'Save Figures & Data');
        set(handles.save_btn, 'BackgroundColor', [0.7, 0.9, 1.0]);
        return;
    end
    
    % Get current timestamp
    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    
    % Save Excel file with calculation results
    if evalin('base', 'exist(''export_data'', ''var'')')
        export_data = evalin('base', 'export_data');
        excel_filename = fullfile(save_dir, [export_data.material_name '_Results_' timestamp '.xlsx']);
        save_results_to_excel(export_data, excel_filename);
        fprintf('Excel file saved: %s\n', excel_filename);
    else
        fprintf('No calculation results available for Excel export.\n');
    end
    
    % Save all figures
    fig_handles = findall(0, 'Type', 'figure');
    for i = 1:length(fig_handles)
        fig_name = sprintf('Figure_%d_%s', get(fig_handles(i), 'Number'), timestamp);
        
        % PNG format (high resolution for publication)
        print(fig_handles(i), fullfile(save_dir, [fig_name '.png']), '-dpng', '-r300');
        
        % EMF format (high quality for publication)
        try
            print(fig_handles(i), fullfile(save_dir, [fig_name '.emf']), '-dmeta', '-r600', '-painters');
        catch
            % EMF not supported on some systems, try alternative
            try
                print(fig_handles(i), fullfile(save_dir, [fig_name '.eps']), '-deps', '-r600', '-painters');
            catch
                % EPS also failed
            end
        end
        
        % FIG format
        saveas(fig_handles(i), fullfile(save_dir, [fig_name '.fig']), 'fig');
    end
    
    % Save workspace
    workspace_file = fullfile(save_dir, ['workspace_' timestamp '.mat']);
    evalin('base', ['save(''' workspace_file ''')']);
    
    % Save the contents of the command window
    command_output_file = fullfile(save_dir, ['command_output_' timestamp '.txt']);
    save_command_window_text(command_output_file);
    
    % Restore the save button state
    set(handles.save_btn, 'Enable', 'on');
    set(handles.save_btn, 'String', 'Save Figures & Data');
    set(handles.save_btn, 'BackgroundColor', [0.7, 0.9, 1.0]);
    

    % Automatically open the saved folder
    try
        if ispc
            % Windows system
            winopen(save_dir);
        elseif ismac
            % macOS system
            system(['open "' save_dir '"']);
        elseif isunix
            % Linux system
            system(['xdg-open "' save_dir '"']);
        end
        fprintf('Automatically opened folder: %s\n', save_dir);
    catch ME
        fprintf('Could not automatically open folder: %s\n', ME.message);
    end

    msgbox(['Results saved to: ' save_dir], 'Save Complete', 'modal');
catch ME
    % Restore button state even when errors occur
    set(handles.save_btn, 'Enable', 'on');
    set(handles.save_btn, 'String', 'Save Figures & Data');
    set(handles.save_btn, 'BackgroundColor', [0.7, 0.9, 1.0]);

    msgbox(['Error saving data: ' ME.message], 'Error', 'error', 'modal');
end

end



function save_command_window_text(filename)
% Save all text in the command window to a file
try
    % Method 1: Get the command window text using the Java interface (applicable to MATLAB 2021b)
    warning('off', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
    
    % Get the Java object of the command window
    jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
    jCmdWin = jDesktop.getClient('Command Window');
    
    if ~isempty(jCmdWin)
        % Get the text component of the command window
        jTextArea = jCmdWin.getComponent(0).getViewport.getComponent(0);
        
        % Get all text
        cmdText = char(jTextArea.getText());
        
        % Writing to a file
        fid = fopen(filename, 'w', 'n', 'UTF-8');
        if fid ~= -1
            fprintf(fid, '%s', cmdText);
            fclose(fid);
            fprintf('Command window text saved to: %s\n', filename);
        else
            error('Cannot create file: %s', filename);
        end
    else
        error('Cannot access command window');
    end
    
catch ME1
    try
        % Method 2: If method 1 fails, use an alternative approach
        fprintf('Method 1 failed, trying alternative method...\n');
        
        % Create an empty file and write the description
        fid = fopen(filename, 'w', 'n', 'UTF-8');
        if fid ~= -1
            fprintf(fid, 'Command Window Text Export\n');
            fprintf(fid, '==========================\n');
            fprintf(fid, 'Export Time: %s\n\n', datestr(now));
            fprintf(fid, 'Note: Automatic text extraction failed.\n');
            fprintf(fid, 'To manually save command window text:\n');
            fprintf(fid, '1. Click in the Command Window\n');
            fprintf(fid, '2. Press Ctrl+A to select all text\n');
            fprintf(fid, '3. Press Ctrl+C to copy\n');
            fprintf(fid, '4. Open this file and paste the content\n\n');
            fprintf(fid, 'Error details: %s\n', ME1.message);
            fclose(fid);
            
            % Display manual operation prompts
            msgbox(['Automatic command window text extraction failed. ' ...
                   'Please manually copy text from Command Window and paste into: ' ...
                   filename], 'Manual Action Required', 'warn');
        end
        
    catch ME2
        fprintf('Both methods failed: %s\n', ME2.message);
        msgbox('Failed to save command window text. Please manually copy and save the text.', 'Error', 'error');
    end
end
end

function close_all_Callback(hObject, eventdata, ~)
% Close all windows EXCEPT the control panel
% Save current parameters before closing figures
fig = ancestor(hObject, 'figure');
handles = get(fig, 'UserData');
save_current_parameters(handles);

all_figs = findall(0, 'Type', 'figure');
control_panel_fig = ancestor(hObject, 'figure');

for i = 1:length(all_figs)
    if all_figs(i) ~= control_panel_fig
        close(all_figs(i));
    end
end
end

function close_panel(hObject, eventdata)
% Custom close function - save parameters before closing
try
    handles = get(hObject, 'UserData');
    if isstruct(handles)
        save_current_parameters(handles);
        fprintf('Parameters saved before closing.\n');
    end

    % Clean up any remaining progress timers
    try
        all_timers = timerfindall('Name', 'ProgressUpdateTimer');
        for i = 1:length(all_timers)
            if isvalid(all_timers(i))
                stop(all_timers(i));
                delete(all_timers(i));
            end
        end
        if ~isempty(all_timers)
            fprintf('Cleaned up %d progress timer(s) before closing.\n', length(all_timers));
        end
    catch
        % Silent cleanup failure
    end
    
    % Clean up parallel pool initialization timer if exists
    try
        pool_init_timers = timerfindall('Name', 'ParallelPoolInit');
        for i = 1:length(pool_init_timers)
            if isvalid(pool_init_timers(i))
                stop(pool_init_timers(i));
                delete(pool_init_timers(i));
            end
        end
    catch
        % Silent cleanup failure
    end
    
    % Close parallel pool when GUI closes
    cleanup_parallel_pool_on_close();

catch
    % Silent save failure
end
delete(hObject);

end

function save_current_parameters(handles)
%% SAVE_CURRENT_PARAMETERS - Persist material-specific settings
%
% PURPOSE:
% Saves current GUI parameter values to file for automatic restoration
% in future sessions, maintaining separate parameter sets for each material.
%
% PARAMETER PERSISTENCE STRATEGY:
% - Material-specific parameter storage prevents cross-contamination
% - Calculation parameters saved per material for workflow efficiency
% - Initial guess coefficients preserved for iterative refinement
% - Session management enables seamless workflow continuation
try
    % Get current material selection
    current_material = get(handles.material_popup, 'Value');
    material_names = {'AZ91', 'Al_4wtCu', 'Al_Cu', 'Al_Fe', 'Sn_Ag', 'Sn_Cu', 'SCN_Acetone', 'Mg_Alloy', 'Custom'};
    material_key = material_names{current_material};
    
    % Load existing material parameters or create new structure
    if exist('material_memory.mat', 'file')
        load('material_memory.mat', 'material_params');
    else
        material_params = struct();
    end
    
    % Create parameter structure for current material
    params = struct();
    
    % Calculation parameters (these are material-specific)
    params.dt_min = get(handles.dt_min_edit, 'String');
    params.dt_max = get(handles.dt_max_edit, 'String');
    params.sampling_interval = get(handles.sampling_interval_edit, 'String');
    params.split_pt = get(handles.split_pt_edit, 'String');
    params.lower_density_mult = get(handles.lower_density_mult_edit, 'String');
    params.upper_density_mult = get(handles.upper_density_mult_edit, 'String');
    params.low_density = get(handles.low_density_edit, 'String');
    params.high_density = get(handles.high_density_edit, 'String');
    params.c0_values = get(handles.c0_values_edit, 'String');
    
    % Initial guess function parameters (always save numerical values)
    params.v_a = get(handles.v_a_edit, 'String');
    params.v_b = get(handles.v_b_edit, 'String');
    params.v_c = get(handles.v_c_edit, 'String');
    params.rd = get(handles.rd_edit, 'String');
    params.re = get(handles.re_edit, 'String');
    params.rf = get(handles.rf_edit, 'String');
    
    % Checkbox settings (material-specific)
    params.include_thermal_undercooling = get(handles.thermal_undercooling_checkbox, 'Value');
    params.use_analytical_approximation = get(handles.use_analytical_checkbox, 'Value');
    
    % Store parameters for this material
    material_params.(material_key) = params;
    
    % Also save the last selected material
    material_params.last_selected_material = current_material;
    
    % Save to file
    save('material_memory.mat', 'material_params');
    fprintf('Parameters saved for material: %s\n', material_key);
    
catch ME
    fprintf('Error saving parameters: %s\n', ME.message);
end
end

function success = load_saved_parameters(handles)
% Load saved parameters for current material
success = false;

try
    if ~exist('material_memory.mat', 'file')
        fprintf('No saved material parameters found, using defaults\n');
        return;
    end
    
    load('material_memory.mat', 'material_params');
    
    % Get last selected material or default to 1
    if isfield(material_params, 'last_selected_material')
        last_material = material_params.last_selected_material;
        set(handles.material_popup, 'Value', last_material);
        fprintf('Restored last selected material: %d\n', last_material);
    else
        last_material = 1;
        set(handles.material_popup, 'Value', last_material);
    end
    
    % Load parameters for the selected material
    load_material_specific_parameters(handles, last_material, material_params);
    
    success = true;
    fprintf('Material-specific parameters loaded successfully\n');
    
catch ME
    fprintf('Error loading parameters: %s\n', ME.message);
    success = false;
end
end

function load_material_specific_parameters(handles, material_index, material_params)
% Load parameters specific to the selected material
material_names = {'AZ91', 'Al_4wtCu', 'Al_Cu', 'Al_Fe', 'Sn_Ag', 'Sn_Cu', 'SCN_Acetone', 'Mg_Alloy', 'Custom'};
material_key = material_names{material_index};

% Check if we have saved parameters for this material
if isfield(material_params, material_key)
    params = material_params.(material_key);
    fprintf('Loading saved parameters for material: %s\n', material_key);
    
    % Load calculation parameters
    if isfield(params, 'dt_min'), set(handles.dt_min_edit, 'String', params.dt_min); end
    if isfield(params, 'dt_max'), set(handles.dt_max_edit, 'String', params.dt_max); end
    if isfield(params, 'sampling_interval'), set(handles.sampling_interval_edit, 'String', params.sampling_interval); end
    if isfield(params, 'split_pt'), set(handles.split_pt_edit, 'String', params.split_pt); end
    if isfield(params, 'lower_density_mult'), set(handles.lower_density_mult_edit, 'String', params.lower_density_mult); end
    if isfield(params, 'upper_density_mult'), set(handles.upper_density_mult_edit, 'String', params.upper_density_mult); end
    if isfield(params, 'low_density'), set(handles.low_density_edit, 'String', params.low_density); end
    if isfield(params, 'high_density'), set(handles.high_density_edit, 'String', params.high_density); end
    if isfield(params, 'c0_values'), set(handles.c0_values_edit, 'String', params.c0_values); end
    
    % Load initial guess parameters
    if isfield(params, 'v_a'), set(handles.v_a_edit, 'String', params.v_a); end
    if isfield(params, 'v_b'), set(handles.v_b_edit, 'String', params.v_b); end
    if isfield(params, 'v_c'), set(handles.v_c_edit, 'String', params.v_c); end
    if isfield(params, 'rd'), set(handles.rd_edit, 'String', params.rd); end
    if isfield(params, 're'), set(handles.re_edit, 'String', params.re); end
    if isfield(params, 'rf'), set(handles.rf_edit, 'String', params.rf); end
    
    % Load checkbox settings
    if isfield(params, 'include_thermal_undercooling')
        set(handles.thermal_undercooling_checkbox, 'Value', params.include_thermal_undercooling);
    end
    if isfield(params, 'use_analytical_approximation')
        set(handles.use_analytical_checkbox, 'Value', params.use_analytical_approximation);
    end
else
    fprintf('No saved parameters for material: %s, using defaults\n', material_key);
    % Will use default material parameters set by material_popup_Callback
end
end






function save_results_to_excel(export_data, filename)
% Save calculation results to Excel file with specified format

try
    if isempty(export_data.all_results)
        fprintf('No calculation results to save.\n');
        return;
    end
    
    % Count valid results
    valid_results = 0;
    max_rows = 0;
    for i = 1:length(export_data.all_results)
        if ~isempty(export_data.all_results{i}) && isfield(export_data.all_results{i}, 'deltaT') && ...
           ~isempty(export_data.all_results{i}.deltaT)
            valid_results = valid_results + 1;
            max_rows = max(max_rows, length(export_data.all_results{i}.deltaT));
        end
    end
    
    if valid_results == 0
        fprintf('No valid calculation results to save.\n');
        return;
    end
    
    fprintf('Preparing Excel export: %d concentrations, maximum %d data points per concentration.\n', ...
            valid_results, max_rows);
    
    % === Sheet 1: Calculation Results ===
    total_cols = valid_results * 3;
    excel_data = cell(max_rows + 2, total_cols);
    
    % Fill headers and data for Sheet 1
    col_idx = 1;
    for i = 1:length(export_data.all_results)
        if ~isempty(export_data.all_results{i}) && isfield(export_data.all_results{i}, 'deltaT') && ...
           ~isempty(export_data.all_results{i}.deltaT)
            
            % First row: Headers
            excel_data{1, col_idx} = 'ΔT (K)';
            excel_data{1, col_idx + 1} = 'R (m)';
            excel_data{1, col_idx + 2} = 'V (m/s)';
            
            % Second row: Concentration values (in second and third columns of each group)
            excel_data{2, col_idx + 1} = export_data.all_results{i}.C0;
            excel_data{2, col_idx + 2} = export_data.all_results{i}.C0;
            
            % Data rows (starting from row 3)
            data_length = length(export_data.all_results{i}.deltaT);
            for row = 1:data_length
                excel_data{row + 2, col_idx} = export_data.all_results{i}.deltaT(row);
                excel_data{row + 2, col_idx + 1} = export_data.all_results{i}.R(row);
                excel_data{row + 2, col_idx + 2} = export_data.all_results{i}.V(row);
            end
            
            col_idx = col_idx + 3;
        end
    end
    
    % Write Sheet 1
    writecell(excel_data, filename, 'Sheet', 'Calculation Results');
    
    % === Sheet 2: Parameters and Fitting Results ===
    sheet2_data = create_parameters_sheet(export_data);
    writecell(sheet2_data, filename, 'Sheet', 'Parameters and Fitting');
    
    fprintf('Excel file successfully saved: %s\n', filename);
    fprintf('File contains %d concentrations with calculation results.\n', valid_results);
    fprintf('Parameters and fitting results saved to Sheet 2.\n');
    
catch ME
    fprintf('Error saving Excel file: %s\n', ME.message);
end

end

function sheet2_data = create_parameters_sheet(export_data)
% Create Sheet 2 with parameters and fitting results

% Initialize cell array for Sheet 2
sheet2_data = {};
row = 1;

% === Material Information ===
sheet2_data{row, 1} = 'MATERIAL INFORMATION';
sheet2_data{row, 1} = format_header_cell(sheet2_data{row, 1});
row = row + 1;

sheet2_data{row, 1} = 'Selected Material:';
sheet2_data{row, 2} = export_data.material_name;
row = row + 2;

% === Material Parameters ===
sheet2_data{row, 1} = 'MATERIAL PARAMETERS';
sheet2_data{row, 1} = format_header_cell(sheet2_data{row, 1});
row = row + 1;

if isfield(export_data, 'material_params')
    mp = export_data.material_params;
    sheet2_data{row, 1} = 'Specific Heat Capacity (Cpv)';
    sheet2_data{row, 2} = sprintf('%.8g J/m³·K', mp.Cpv);
    row = row + 1;
    
    sheet2_data{row, 1} = 'Latent Heat (DHv)';
    sheet2_data{row, 2} = sprintf('%.8g J/m³', mp.DHv);
    row = row + 1;
    
    sheet2_data{row, 1} = 'Thermal Diffusivity (a)';
    sheet2_data{row, 2} = sprintf('%.8g m²/s', mp.a);
    row = row + 1;
    
    sheet2_data{row, 1} = 'Gibbs-Thomson Coefficient';
    sheet2_data{row, 2} = sprintf('%.8g K·m', mp.Gibbs_Tom);
    row = row + 1;
    
    sheet2_data{row, 1} = 'Diffusion Coefficient (D_Al)';
    sheet2_data{row, 2} = sprintf('%.8g m²/s', mp.D_Al);
    row = row + 1;
    
    sheet2_data{row, 1} = 'Liquidus Slope (m_Al)';
    sheet2_data{row, 2} = sprintf('%.8g K/wt%%', mp.m_Al);
    row = row + 1;
    
    sheet2_data{row, 1} = 'Partition Coefficient (k0_Al)';
    sheet2_data{row, 2} = sprintf('%.8g', mp.k0_Al);
    row = row + 1;
    
    sheet2_data{row, 1} = 'Stability Constant (σ)';
    sheet2_data{row, 2} = sprintf('%.8g', mp.sigma);
    row = row + 1;
end
row = row + 1;

% === Calculation Parameters ===
sheet2_data{row, 1} = 'CALCULATION PARAMETERS';
sheet2_data{row, 1} = format_header_cell(sheet2_data{row, 1});
row = row + 1;

if isfield(export_data, 'calc_params')
    cp = export_data.calc_params;
    sheet2_data{row, 1} = 'ΔT minimum';
    sheet2_data{row, 2} = sprintf('%.8g K', cp.deltaT_min);
    row = row + 1;
    
    sheet2_data{row, 1} = 'ΔT maximum';
    sheet2_data{row, 2} = sprintf('%.8g K', cp.deltaT_max);
    row = row + 1;
    
    sheet2_data{row, 1} = 'Sampling Interval';
    sheet2_data{row, 2} = sprintf('%.8g', cp.sampling_interval);
    row = row + 1;
    
    sheet2_data{row, 1} = 'Split Point';
    sheet2_data{row, 2} = sprintf('%.8g', cp.split_point);
    row = row + 1;
    
    sheet2_data{row, 1} = 'Lower Density Multiplier';
    sheet2_data{row, 2} = sprintf('%.8g', cp.lower_density_mult);
    row = row + 1;
    
    sheet2_data{row, 1} = 'Upper Density Multiplier';
    sheet2_data{row, 2} = sprintf('%.8g', cp.upper_density_mult);
    row = row + 1;
    
    sheet2_data{row, 1} = 'Low End Density';
    sheet2_data{row, 2} = sprintf('%.8g', cp.low_end_density);
    row = row + 1;
    
    sheet2_data{row, 1} = 'High End Density';
    sheet2_data{row, 2} = sprintf('%.8g', cp.high_end_density);
    row = row + 1;
    
    sheet2_data{row, 1} = 'C0 Values';
    sheet2_data{row, 2} = cp.c0_values;
    row = row + 1;
end
row = row + 1;

% === Calculation Time ===
sheet2_data{row, 1} = 'CALCULATION TIME';
sheet2_data{row, 1} = format_header_cell(sheet2_data{row, 1});
row = row + 1;

if isfield(export_data, 'calculation_time')
    sheet2_data{row, 1} = 'Total Computation Time';
    sheet2_data{row, 2} = sprintf('%.2f seconds', export_data.calculation_time);
    row = row + 1;
end
row = row + 1;

% === Initial Guess Functions ===
sheet2_data{row, 1} = 'INITIAL GUESS FUNCTIONS';
sheet2_data{row, 1} = format_header_cell(sheet2_data{row, 1});
row = row + 1;

if isfield(export_data, 'calc_params')
    cp = export_data.calc_params;
    sheet2_data{row, 1} = 'V_guess Function';
    sheet2_data{row, 2} = sprintf('V = (ΔT/(%.8g × C0^%.8g))^%.8g', cp.V_a, cp.V_b, cp.V_c);
    row = row + 1;
    
    sheet2_data{row, 1} = 'V_a';
    sheet2_data{row, 2} = sprintf('%.8g', cp.V_a);
    row = row + 1;
    
    sheet2_data{row, 1} = 'V_b';
    sheet2_data{row, 2} = sprintf('%.8g', cp.V_b);
    row = row + 1;
    
    sheet2_data{row, 1} = 'V_c';
    sheet2_data{row, 2} = sprintf('%.8g', cp.V_c);
    row = row + 1;
    
    sheet2_data{row, 1} = 'R_guess Function';
    sheet2_data{row, 2} = sprintf('R = %.8g × C0^%.8g × ΔT^%.8g', cp.R_d, cp.R_e, cp.R_f);
    row = row + 1;
    
    sheet2_data{row, 1} = 'R_d';
    sheet2_data{row, 2} = sprintf('%.8g', cp.R_d);
    row = row + 1;
    
    sheet2_data{row, 1} = 'R_e';
    sheet2_data{row, 2} = sprintf('%.8g', cp.R_e);
    row = row + 1;
    
    sheet2_data{row, 1} = 'R_f';
    sheet2_data{row, 2} = sprintf('%.8g', cp.R_f);
    row = row + 1;
end
row = row + 1;

% === Surface Fitting Results from Figures 8, 9, 10 ===
sheet2_data{row, 1} = 'SURFACE FITTING RESULTS';
sheet2_data{row, 1} = format_header_cell(sheet2_data{row, 1});
row = row + 1;

% Figure 8: V Surface Fitting
try
    if evalin('base', 'exist(''V_surface_fitting'', ''var'')')
        V_fit = evalin('base', 'V_surface_fitting');
        sheet2_data{row, 1} = sprintf('Figure %d - V Surface Fitting', V_fit.figure_number);
        row = row + 1;
        
        sheet2_data{row, 1} = 'Equation';
        sheet2_data{row, 2} = V_fit.equation;
        row = row + 1;
        
        sheet2_data{row, 1} = 'Complete Expression';
        sheet2_data{row, 2} = sprintf('V = (ΔT/(%.8g×C0^%.8g))^%.8g', V_fit.a, V_fit.b, V_fit.c);
        row = row + 1;
        
        sheet2_data{row, 1} = 'Parameter a';
        sheet2_data{row, 2} = sprintf('%.8g', V_fit.a);
        row = row + 1;
        
        sheet2_data{row, 1} = 'Parameter b';
        sheet2_data{row, 2} = sprintf('%.8g', V_fit.b);
        row = row + 1;
        
        sheet2_data{row, 1} = 'Parameter c';
        sheet2_data{row, 2} = sprintf('%.8g', V_fit.c);
        row = row + 1;
        
        sheet2_data{row, 1} = 'R-squared';
        sheet2_data{row, 2} = sprintf('%.8g', V_fit.R_squared);
        row = row + 1;
        row = row + 1;
    else
        sheet2_data{row, 1} = 'Figure 8 - V Surface Fitting: Not Available';
        row = row + 2;
    end
catch ME
    sheet2_data{row, 1} = 'Figure 8 - V Surface Fitting: Error';
    sheet2_data{row, 2} = ME.message;
    row = row + 2;
end

% Figure 9: R Surface Fitting
try
    if evalin('base', 'exist(''R_surface_fitting'', ''var'')')
        R_fit = evalin('base', 'R_surface_fitting');
        sheet2_data{row, 1} = sprintf('Figure %d - R Surface Fitting', R_fit.figure_number);
        row = row + 1;
        
        sheet2_data{row, 1} = 'Equation';
        sheet2_data{row, 2} = R_fit.equation;
        row = row + 1;
        
        sheet2_data{row, 1} = 'Complete Expression';
        sheet2_data{row, 2} = sprintf('R = %.8g×C0^%.8g×ΔT^%.8g', R_fit.d, R_fit.e, R_fit.f);
        row = row + 1;
        
        sheet2_data{row, 1} = 'Parameter d';
        sheet2_data{row, 2} = sprintf('%.8g', R_fit.d);
        row = row + 1;
        
        sheet2_data{row, 1} = 'Parameter e';
        sheet2_data{row, 2} = sprintf('%.8g', R_fit.e);
        row = row + 1;
        
        sheet2_data{row, 1} = 'Parameter f';
        sheet2_data{row, 2} = sprintf('%.8g', R_fit.f);
        row = row + 1;
        
        sheet2_data{row, 1} = 'R-squared';
        sheet2_data{row, 2} = sprintf('%.8g', R_fit.R_squared);
        row = row + 1;
        row = row + 1;
    else
        sheet2_data{row, 1} = 'Figure 9 - R Surface Fitting: Not Available';
        row = row + 2;
    end
catch ME
    sheet2_data{row, 1} = 'Figure 9 - R Surface Fitting: Error';
    sheet2_data{row, 2} = ME.message;
    row = row + 2;
end

% Figure 10: ΔT Surface Fitting
try
    if evalin('base', 'exist(''deltaT_surface_fitting'', ''var'')')
        deltaT_fit = evalin('base', 'deltaT_surface_fitting');
        sheet2_data{row, 1} = sprintf('Figure %d - ΔT Surface Fitting', deltaT_fit.figure_number);
        row = row + 1;
        
        sheet2_data{row, 1} = 'Equation';
        sheet2_data{row, 2} = deltaT_fit.equation;
        row = row + 1;
        
        sheet2_data{row, 1} = 'Complete Expression';
        sheet2_data{row, 2} = sprintf('ΔT = %.8g×C0^%.8g×V^%.8g', deltaT_fit.a, deltaT_fit.b, deltaT_fit.c);
        row = row + 1;
        
        sheet2_data{row, 1} = 'Parameter a';
        sheet2_data{row, 2} = sprintf('%.8g', deltaT_fit.a);
        row = row + 1;
        
        sheet2_data{row, 1} = 'Parameter b';
        sheet2_data{row, 2} = sprintf('%.8g', deltaT_fit.b);
        row = row + 1;
        
        sheet2_data{row, 1} = 'Parameter c';
        sheet2_data{row, 2} = sprintf('%.8g', deltaT_fit.c);
        row = row + 1;
        
        sheet2_data{row, 1} = 'R-squared';
        sheet2_data{row, 2} = sprintf('%.8g', deltaT_fit.R_squared);
        row = row + 1;
    else
        sheet2_data{row, 1} = 'Figure 10 - ΔT Surface Fitting: Not Available';
        row = row + 1;
    end
catch ME
    sheet2_data{row, 1} = 'Figure 10 - ΔT Surface Fitting: Error';
    sheet2_data{row, 2} = ME.message;
    row = row + 1;
end

end

function header_text = format_header_cell(text)
% Format header text (this is a placeholder - Excel formatting is limited in basic writecell)
header_text = ['=== ' text ' ==='];
end

function auto_guess_Callback(hObject, eventdata, ~)
% Smart search based on current parameters, with progress shown on the control panel
fig = ancestor(hObject, 'figure');
handles = get(fig, 'UserData');

try
    % Completely suppress output from calculate_VR

    global SUPPRESS_OUTPUT;
    original_suppress_state = SUPPRESS_OUTPUT;
    SUPPRESS_OUTPUT = true;
    
    % Initialise the search status
    set(handles.auto_guess_btn, 'Enable', 'off');
    set(handles.auto_guess_btn, 'String', 'Initializing...');
    set(handles.auto_guess_btn, 'BackgroundColor', [0.9, 0.9, 0.9]);
    
    % Clear the progress display
    reset_progress_display(handles);
    drawnow;
    
    % Retrieve material parameters and test conditions
    material_params = get_material_params_from_gui(handles);
    c0_str = get(handles.c0_values_edit, 'String');
    C0_test_values = eval(c0_str);
    C0_test = C0_test_values(ceil(length(C0_test_values)/2));
    
    % Generate a test deltaT sequence
    deltaT_sequence = get_test_deltaT_sequence(handles);
    deltaT_test = deltaT_sequence(1:min(3, length(deltaT_sequence)));
    
    % Create a search range based on current parameters
    search_ranges = create_search_ranges_from_current(handles);
    total_combinations = calculate_total_combinations(search_ranges);
    
    % Display search information in the command line
    fprintf('\n=== Parameter Search Started ===\n');
    fprintf('Material: %s\n', get_material_name(handles));
    fprintf('Current parameters as center: V_a=%.3f, V_b=%.3f, V_c=%.3f, R_d=%.2e, R_e=%.3f, R_f=%.3f\n', ...
            str2double(get(handles.v_a_edit, 'String')), str2double(get(handles.v_b_edit, 'String')), ...
            str2double(get(handles.v_c_edit, 'String')), str2double(get(handles.rd_edit, 'String')), ...
            str2double(get(handles.re_edit, 'String')), str2double(get(handles.rf_edit, 'String')));
    fprintf('Total combinations: %d\n', total_combinations);
    
    % Check the availability of parallel computing
    use_parallel = setup_parallel_if_available();
    if use_parallel
        fprintf('Using parallel computing (%d workers)\n', gcp().NumWorkers);
    else
        fprintf('Using serial computing\n');
    end
    
    % Execute the search
    if use_parallel
        best_params = parallel_parameter_search_with_progress(search_ranges, C0_test, deltaT_test, ...
                                                            material_params, handles, total_combinations);
    else
        best_params = serial_parameter_search_with_progress(search_ranges, C0_test, deltaT_test, ...
                                                          material_params, handles, total_combinations);
    end
    
    % Apply the results
    if best_params.score > 0
        apply_best_parameters(handles, best_params);
        fprintf('=== Search Completed Successfully ===\n');
        fprintf('Best score: %.2f (out of %.0f test points)\n', best_params.score, length(deltaT_test));
        
        success_msg = sprintf('Search completed!\nTested %d combinations\nBest score: %.2f/%d\nParameters updated.', ...
                            total_combinations, best_params.score, length(deltaT_test));
        msgbox(success_msg, 'Success');
    else
        fprintf('=== Search Failed ===\n');
        msgbox('No better parameters found. Current parameters may already be optimal.', 'Info');
    end
    
catch ME
    fprintf('=== Search Error ===\n');
    fprintf('Error: %s\n', ME.message);
    msgbox(['Error during search: ' ME.message], 'Error', 'error');
end

% Restore the previous state
SUPPRESS_OUTPUT = original_suppress_state;
set(handles.auto_guess_btn, 'Enable', 'on');
set(handles.auto_guess_btn, 'String', 'Auto Find Parameters');
set(handles.auto_guess_btn, 'BackgroundColor', [0.8, 1.0, 0.8]);
end



function apply_best_parameters(handles, best_params)
% Apply the best parameters to the GUI
set(handles.v_a_edit, 'String', num2str(best_params.V_a));
set(handles.v_b_edit, 'String', num2str(best_params.V_b));
set(handles.v_c_edit, 'String', num2str(best_params.V_c));
set(handles.rd_edit, 'String', num2str(best_params.R_d));
set(handles.re_edit, 'String', num2str(best_params.R_e));
set(handles.rf_edit, 'String', num2str(best_params.R_f));
end

function material_params = get_material_params_from_gui(handles)
% Extract material parameters from the GUI
material_params = struct();
material_params.Cpv = str2double(get(handles.cp_edit, 'String'));
material_params.DHv = str2double(get(handles.dh_edit, 'String'));
material_params.a = str2double(get(handles.a_edit, 'String'));
material_params.Gibbs_Tom = str2double(get(handles.gibbs_edit, 'String'));
material_params.D_Al = str2double(get(handles.dal_edit, 'String'));
material_params.m_Al = str2double(get(handles.mal_edit, 'String'));
material_params.k0_Al = str2double(get(handles.k0al_edit, 'String'));
material_params.sigma = str2double(get(handles.sigma_edit, 'String'));
end

function score = evaluate_parameter_combination(params, C0_test, deltaT_test, material_params)
% Evaluate a single parameter combination with all output suppressed
V_a = params(1); V_b = params(2); V_c = params(3);
R_d = params(4); R_e = params(5); R_f = params(6);

try
    % Create the guess functions
    V_guess_func = @(deltaT) (deltaT./(V_a * C0_test.^V_b)).^V_c;
    R_guess_func = @(deltaT) R_d * (C0_test.^R_e) .* (deltaT.^R_f);
    
    % Temporarily redirect output to suppress command-line display
    orig_state = warning('off', 'all');
    
    % Call the computation function
    [result, ~, ~] = calculate_VR(V_guess_func, R_guess_func, C0_test, deltaT_test, material_params);
    
    % Restore the warning state
    warning(orig_state);
    
    % Calculate the score
    score = size(result, 1) / length(deltaT_test);
    
catch
    score = 0;
end
end


function material_name = get_material_name(handles)
% Get the name of the currently selected material
material_names = {'AZ91', 'Al-4wt%Cu', 'Al-Cu', 'Al-Fe', 'Sn-Ag', 'Sn-Cu', 'SCN-Acetone(LGK 1984)', 'Mg-Alloy(Lin 2009)', 'Custom'};
material_selection = get(handles.material_popup, 'Value');
material_name = material_names{material_selection};
end


function score = evaluate_parameter_combination_quiet(params, C0_test, deltaT_test, material_params)
% Completely silent parameter evaluation function
V_a = params(1); V_b = params(2); V_c = params(3);
R_d = params(4); R_e = params(5); R_f = params(6);

try
    % Create the guess functions
    V_guess_func = @(deltaT) (deltaT./(V_a * C0_test.^V_b)).^V_c;
    R_guess_func = @(deltaT) R_d * (C0_test.^R_e) .* (deltaT.^R_f);
    
    % Temporarily suppress all output and warnings
    orig_warning_state = warning('off', 'all');
    
    % Redirect fprintf output to a null device (completely silent)
    global SUPPRESS_OUTPUT;
    original_suppress = SUPPRESS_OUTPUT;
    SUPPRESS_OUTPUT = true;
    
    % Calling the calculation function
    [result, ~, ~] = calculate_VR(V_guess_func, R_guess_func, C0_test, deltaT_test, material_params);
    
    % Restore the previous state
    SUPPRESS_OUTPUT = original_suppress;
    warning(orig_warning_state);
    
    % Calculate the score
    score = size(result, 1) / length(deltaT_test);
    
catch
    score = 0;
end
end


function search_ranges = create_search_ranges_from_current(handles)
% Create a search range based on current parameters (reduce points for speed)
current_V_a = str2double(get(handles.v_a_edit, 'String'));
current_V_b = str2double(get(handles.v_b_edit, 'String'));
current_V_c = str2double(get(handles.v_c_edit, 'String'));
current_R_d = str2double(get(handles.rd_edit, 'String'));
current_R_e = str2double(get(handles.re_edit, 'String'));
current_R_f = str2double(get(handles.rf_edit, 'String'));

% Reduce the number of search points (3×3×3×3×3×3 = 729 combinations, faster)
search_ranges.V_a = linspace(current_V_a * 0.7, current_V_a * 1.3, 3);
search_ranges.V_b = linspace(max(0.1, current_V_b * 0.7), current_V_b * 1.3, 3);
search_ranges.V_c = linspace(max(1.0, current_V_c * 0.7), current_V_c * 1.3, 3);
search_ranges.R_d = linspace(current_R_d * 0.2, current_R_d * 5, 3);
search_ranges.R_e = linspace(current_R_e - 0.5, current_R_e + 0.5, 3);
search_ranges.R_f = linspace(current_R_f - 0.5, current_R_f + 0.5, 3);
end

function deltaT_sequence = get_test_deltaT_sequence(handles)
% Obtain the deltaT sequence for testing
deltaT_min = str2double(get(handles.dt_min_edit, 'String'));
deltaT_max = str2double(get(handles.dt_max_edit, 'String'));
sampling_interval = str2double(get(handles.sampling_interval_edit, 'String'));
split_point = str2double(get(handles.split_pt_edit, 'String'));
lower_density_mult = str2double(get(handles.lower_density_mult_edit, 'String'));
upper_density_mult = str2double(get(handles.upper_density_mult_edit, 'String'));
low_end_density = str2double(get(handles.low_density_edit, 'String'));
high_end_density = str2double(get(handles.high_density_edit, 'String'));

deltaT_sequence = create_deltaT_sequence('deltaT_min', deltaT_min, 'deltaT_max', deltaT_max, ...
                                       'sampling_interval', sampling_interval, 'split_point', split_point, ...
                                       'lower_density_mult', lower_density_mult, 'upper_density_mult', upper_density_mult, ...
                                       'low_end_density', low_end_density, 'high_end_density', high_end_density);
end

function total = calculate_total_combinations(search_ranges)
% Calculate the total number of combinations
total = length(search_ranges.V_a) * length(search_ranges.V_b) * length(search_ranges.V_c) * ...
        length(search_ranges.R_d) * length(search_ranges.R_e) * length(search_ranges.R_f);
end

function use_parallel = setup_parallel_if_available()
% Set up parallel computing (if available)
use_parallel = false;
try
    if license('test', 'Distrib_Computing_Toolbox')
        if isempty(gcp('nocreate'))
            % Automatically detect the number of cores, use total cores minus two, with a minimum of two cores
            total_cores = feature('numcores');
            worker_count = max(2, total_cores - 2);
            parpool('local', worker_count);
            fprintf('Detected %d CPU cores, using %d workers for parallel computing\n', total_cores, worker_count);
        end
        use_parallel = true;
    end
catch
    % Parallel computing not available
end
end

function reset_progress_display(handles)
% Reset the progress display
set(handles.progress_text_bar, 'String', '[░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░]');
set(handles.progress_percent, 'String', '0%');
set(handles.progress_elapsed, 'String', 'Elapsed: 0.00 min');
set(handles.progress_remaining, 'String', 'Remaining: -- min');
end

function update_progress_display_simple(handles, completed, total, elapsed_time, best_score)
% Simple update of the progress display
try
    if ~isvalid(handles.progress_text_bar)
        return;
    end
    
    % Update the progress bar
    progress_percent = (completed / total) * 100;
    filled_chars = round((completed / total) * 40);
    empty_chars = 40 - filled_chars;
    progress_bar = ['[' repmat('█', 1, filled_chars) repmat('░', 1, empty_chars) ']'];
    
    set(handles.progress_text_bar, 'String', progress_bar);
    set(handles.progress_percent, 'String', sprintf('%.1f%%', progress_percent));
    set(handles.progress_elapsed, 'String', sprintf('Elapsed: %.2f min', elapsed_time));
    
    % Estimate the remaining time
    if completed > 0 && elapsed_time > 0
        remaining_time = (elapsed_time * total / completed) - elapsed_time;
        set(handles.progress_remaining, 'String', sprintf('Remaining: %.2f min', remaining_time));
    end
    
    drawnow limitrate;
catch
    % Silently handle any errors
end
end

function best_params = serial_parameter_search_with_progress(search_ranges, C0_test, deltaT_test, material_params, handles, total_combinations)
% Serial search with progress display
best_params = struct('score', 0);
current_combination = 0;
start_time = tic;

for V_a = search_ranges.V_a
    for V_b = search_ranges.V_b
        for V_c = search_ranges.V_c
            for R_d = search_ranges.R_d
                for R_e = search_ranges.R_e
                    for R_f = search_ranges.R_f
                        current_combination = current_combination + 1;
                        
                        % Update the progress display
                        if mod(current_combination, max(1, floor(total_combinations/100))) == 0 || current_combination == total_combinations
                            elapsed_time = toc(start_time) / 60;
                            update_progress_display_simple(handles, current_combination, total_combinations, elapsed_time, best_params.score);
                        end
                        
                        % Evaluate the parameters
                        score = evaluate_parameter_combination_quiet([V_a, V_b, V_c, R_d, R_e, R_f], ...
                                                                   C0_test, deltaT_test, material_params);
                        
                        if score > best_params.score
                            best_params.V_a = V_a; best_params.V_b = V_b; best_params.V_c = V_c;
                            best_params.R_d = R_d; best_params.R_e = R_e; best_params.R_f = R_f;
                            best_params.score = score;
                        end
                    end
                end
            end
        end
    end
end

% Final update
elapsed_time = toc(start_time) / 60;
update_progress_display_simple(handles, total_combinations, total_combinations, elapsed_time, best_params.score);
end

function best_params = parallel_parameter_search_with_progress(search_ranges, C0_test, deltaT_test, material_params, handles, total_combinations)
% Parallel search with progress display
best_params = struct('score', 0);

% Create the parameter combination matrix
param_combinations = [];
for V_a = search_ranges.V_a
    for V_b = search_ranges.V_b
        for V_c = search_ranges.V_c
            for R_d = search_ranges.R_d
                for R_e = search_ranges.R_e
                    for R_f = search_ranges.R_f
                        param_combinations = [param_combinations; V_a, V_b, V_c, R_d, R_e, R_f];
                    end
                end
            end
        end
    end
end

% Process in parallel batches
batch_size = 100;
num_batches = ceil(size(param_combinations, 1) / batch_size);
start_time = tic;

for batch = 1:num_batches
    start_idx = (batch-1) * batch_size + 1;
    end_idx = min(batch * batch_size, size(param_combinations, 1));
    current_batch = param_combinations(start_idx:end_idx, :);
    
    % Evaluate the current batch in parallel
    batch_scores = zeros(size(current_batch, 1), 1);
    parfor i = 1:size(current_batch, 1)
        batch_scores(i) = evaluate_parameter_combination_quiet(current_batch(i, :), C0_test, deltaT_test, material_params);
    end
    
    % Update the best parameters
    [batch_best_score, batch_best_idx] = max(batch_scores);
    if batch_best_score > best_params.score
        best_combo = current_batch(batch_best_idx, :);
        best_params.V_a = best_combo(1); best_params.V_b = best_combo(2); best_params.V_c = best_combo(3);
        best_params.R_d = best_combo(4); best_params.R_e = best_combo(5); best_params.R_f = best_combo(6);
        best_params.score = batch_best_score;
    end
    
    % Update the progress
    elapsed_time = toc(start_time) / 60;
    update_progress_display_simple(handles, end_idx, total_combinations, elapsed_time, best_params.score);
end
end

function stop_calculation_Callback(hObject, eventdata, ~)
% Stop calculation and clean up parallel pool
global STOP_CALCULATION;

try
    % Set stop flag
    STOP_CALCULATION = true;
    
    % Get handles
    fig = ancestor(hObject, 'figure');
    handles = get(fig, 'UserData');
    
    % Disable stop button
    set(handles.stop_btn, 'Enable', 'off');
    set(handles.stop_btn, 'String', 'Stopping...');
    
    fprintf('Stop signal sent. Cleaning up parallel computation...\n');
    
    % Clean up parallel pool
    cleanup_parallel_pool();
    
    % Re-enable start button
    set(handles.start_btn, 'Enable', 'on');
    set(handles.start_btn, 'String', 'Start Calculation');
    set(handles.start_btn, 'BackgroundColor', [0.7, 1.0, 0.7]);
    
    % Reset stop button
    set(handles.stop_btn, 'String', 'Stop Calculation');
    
    % Hide instruction text
    if isfield(handles, 'instruction_text')
        set(handles.instruction_text, 'Visible', 'off');
    end
    
    msgbox('Calculation stopped and parallel pool cleaned up.', 'Stopped', 'modal');
    
catch ME
    fprintf('Error during stop: %s\n', ME.message);
end
end

function [V_guess, R_guess] = get_SCN_fitted_functions(C0, deltaT)
% GET_SCN_FITTED_FUNCTIONS - SCN-Acetone interpolation functions based on sixth-order polynomial fitting
% Use sixth-order polynomials to fit V and R

% Polynomial coefficients for V - sixth-order (from highest to lowest degree)
V_05K_poly_coeffs = [-6.659290507617004e-04, 2.369698931357119e-03, -3.396311141455152e-03, 2.487586100839964e-03, -9.373914269260196e-04, 1.173158416800813e-04, 3.764308570717425e-05];
V_09K_poly_coeffs = [-5.533479076924480e-03, 1.958381856164818e-02, -2.749709894968078e-02, 1.941241773272147e-02, -7.029175352673138e-03, 9.566442665513241e-04, 1.691807667749731e-04];

% Polynomial coefficients for R - sixth-order (from highest to lowest degree)
R_05K_poly_coeffs = [3.933800623029993e-04, -1.399829155419792e-03, 1.993583146697738e-03, -1.455094212818860e-03, 5.807424879235848e-04, -1.226741862485822e-04, 2.463444587967703e-05];
R_09K_poly_coeffs = [4.062364954450699e-04, -1.402966115638020e-03, 1.906008843886529e-03, -1.294100697718033e-03, 4.630895797590333e-04, -8.404160371893581e-05, 1.193829099220079e-05];

% Calculate the values of V and R
V_at_05K = polyval(V_05K_poly_coeffs, C0);
V_at_09K = polyval(V_09K_poly_coeffs, C0);
R_at_05K = polyval(R_05K_poly_coeffs, C0);
R_at_09K = polyval(R_09K_poly_coeffs, C0);

% Linear interpolation
deltaT_clamped = max(0.5, min(0.9, deltaT));
weight_05K = (0.9 - deltaT_clamped) / 0.4;
weight_09K = (deltaT_clamped - 0.5) / 0.4;

V_guess = weight_05K * V_at_05K + weight_09K * V_at_09K;
R_guess = weight_05K * R_at_05K + weight_09K * R_at_09K;

% Ensure positive values
V_guess = max(V_guess, 1e-10);
R_guess = max(R_guess, 1e-10);

end

function thermal_undercooling_Callback(hObject, eventdata, ~)
% Callback for thermal undercooling checkbox
fig = ancestor(hObject, 'figure');
handles = get(fig, 'UserData');

% Update figure's UserData
set(fig, 'UserData', handles);
end

function split_point_edit_Callback(hObject, eventdata, handles)
%% SPLIT_POINT_EDIT_CALLBACK - Validate deltaT sequence split point
%
% PURPOSE:
% Validates split point for two-region deltaT sequence generation and
% enables/disables density control parameters accordingly.
%
% VALIDATION LOGIC:
% Split point must satisfy: deltaT_min < split_point < deltaT_max
% - Valid split point: Enables density multiplier controls for two-region mode
% - Invalid split point: Disables density controls, uses single linear sequence
%
% TWO-REGION SEQUENCE BENEFITS:
% - Higher density sampling in critical undercooling ranges
% - Computational efficiency through adaptive point distribution
% - Enhanced resolution where dendrite behaviour changes rapidly
if nargin < 3
    fig = ancestor(hObject, 'figure');
    handles = get(fig, 'UserData');
end

% Get current values
deltaT_min = str2double(get(handles.dt_min_edit, 'String'));
deltaT_max = str2double(get(handles.dt_max_edit, 'String'));
split_point = str2double(get(handles.split_pt_edit, 'String'));

% Check if values are valid numbers
if isnan(deltaT_min) || isnan(deltaT_max) || isnan(split_point)
    return; % Invalid input, don't change interface
end

% Check if split point is within range (exclusive)
if split_point > deltaT_min && split_point < deltaT_max
    % Split point is valid - enable density controls
    set(handles.lower_density_mult_edit, 'Enable', 'on');
    set(handles.upper_density_mult_edit, 'Enable', 'on');
    set(handles.low_density_edit, 'Enable', 'on');
    set(handles.high_density_edit, 'Enable', 'on');
    
    % Change background color to indicate enabled state
    set(handles.lower_density_mult_edit, 'BackgroundColor', 'white');
    set(handles.upper_density_mult_edit, 'BackgroundColor', 'white');
    set(handles.low_density_edit, 'BackgroundColor', 'white');
    set(handles.high_density_edit, 'BackgroundColor', 'white');
    
    fprintf('Split point %.2f is valid (%.2f < %.2f < %.2f) - Density controls enabled\n', ...
            split_point, deltaT_min, split_point, deltaT_max);
else
    % Split point is invalid - disable density controls
    set(handles.lower_density_mult_edit, 'Enable', 'off');
    set(handles.upper_density_mult_edit, 'Enable', 'off');
    set(handles.low_density_edit, 'Enable', 'off');
    set(handles.high_density_edit, 'Enable', 'off');
    
    % Change background color to indicate disabled state
    set(handles.lower_density_mult_edit, 'BackgroundColor', [0.9, 0.9, 0.9]);
    set(handles.upper_density_mult_edit, 'BackgroundColor', [0.9, 0.9, 0.9]);
    set(handles.low_density_edit, 'BackgroundColor', [0.9, 0.9, 0.9]);
    set(handles.high_density_edit, 'BackgroundColor', [0.9, 0.9, 0.9]);
    
    fprintf('Split point %.2f is invalid (not in range %.2f to %.2f) - Density controls disabled\n', ...
            split_point, deltaT_min, deltaT_max);
end

end

function deltaT_range_edit_Callback(hObject, eventdata, handles)
% Callback for deltaT_min or deltaT_max changes - revalidate split point
if nargin < 3
    fig = ancestor(hObject, 'figure');
    handles = get(fig, 'UserData');
end

% Call split point validation
split_point_edit_Callback(handles.split_pt_edit, [], handles);
end

function analytical_checkbox_Callback(hObject, eventdata, ~)
%% ANALYTICAL_CHECKBOX_CALLBACK - Toggle initial guess function mode
%
% PURPOSE:
% Switches between parametric guess functions and analytical approximations
% based on equations (8.91) and (8.92) from solidification theory literature.
%
% APPROXIMATION MODES:
% - Parametric: V = (ΔT/(a×C₀^b))^c, R = d×C₀^e×ΔT^f
% - Analytical: Uses theoretical expressions from Kurz & Fisher equations
%   V ∝ (D_l×ΔT^2.5)/(Γ×C₀^1.5×(-m×(1-k₀))^1.5)
%   R ∝ (Γ×C₀^0.25×(-m×(1-k₀))^0.25)/ΔT^1.25
%
% CONTROL LOGIC:
% - Disables parameter editing when analytical mode active
% - Maintains material-specific behaviour for SCN-Acetone
% - Provides fallback to parametric mode if analytical functions unavailable

fig = ancestor(hObject, 'figure');
handles = get(fig, 'UserData');

use_analytical = get(hObject, 'Value');
current_material = get(handles.material_popup, 'Value');

if use_analytical
    % Disable parameter edit boxes when using analytical approximation
    set(handles.v_a_edit, 'Enable', 'off');
    set(handles.v_b_edit, 'Enable', 'off');
    set(handles.v_c_edit, 'Enable', 'off');
    set(handles.rd_edit, 'Enable', 'off');
    set(handles.re_edit, 'Enable', 'off');
    set(handles.rf_edit, 'Enable', 'off');
    
    % Change background color to indicate disabled state
    set(handles.v_a_edit, 'BackgroundColor', [0.9, 0.9, 0.9]);
    set(handles.v_b_edit, 'BackgroundColor', [0.9, 0.9, 0.9]);
    set(handles.v_c_edit, 'BackgroundColor', [0.9, 0.9, 0.9]);
    set(handles.rd_edit, 'BackgroundColor', [0.9, 0.9, 0.9]);
    set(handles.re_edit, 'BackgroundColor', [0.9, 0.9, 0.9]);
    set(handles.rf_edit, 'BackgroundColor', [0.9, 0.9, 0.9]);
    
    fprintf('Analytical approximation enabled: Using Eqs. (8.91) & (8.92)\n');
else
    % Re-enable parameter edit boxes (unless SCN material)
    if current_material ~= 7  % Not SCN-Acetone
        set(handles.v_a_edit, 'Enable', 'on');
        set(handles.v_b_edit, 'Enable', 'on');
        set(handles.v_c_edit, 'Enable', 'on');
        set(handles.rd_edit, 'Enable', 'on');
        set(handles.re_edit, 'Enable', 'on');
        set(handles.rf_edit, 'Enable', 'on');
        
        % Restore background color
        set(handles.v_a_edit, 'BackgroundColor', 'white');
        set(handles.v_b_edit, 'BackgroundColor', 'white');
        set(handles.v_c_edit, 'BackgroundColor', 'white');
        set(handles.rd_edit, 'BackgroundColor', 'white');
        set(handles.re_edit, 'BackgroundColor', 'white');
        set(handles.rf_edit, 'BackgroundColor', 'white');
    end
    
    fprintf('Analytical approximation disabled: Using parametric functions\n');
end

% Update figure's UserData
set(fig, 'UserData', handles);

% Do NOT save parameters here - only save when explicitly requested
end

function enable_all_controls(handles)
% Enable all controls for custom material
control_fields = {'cp_edit', 'dh_edit', 'a_edit', 'gibbs_edit', 'dal_edit', 'mal_edit', 'k0al_edit', 'sigma_edit', ...
                  'dt_min_edit', 'dt_max_edit', 'sampling_interval_edit', 'split_pt_edit', ...
                  'lower_density_mult_edit', 'upper_density_mult_edit', 'low_density_edit', 'high_density_edit', ...
                  'c0_values_edit', 'v_a_edit', 'v_b_edit', 'v_c_edit', 'rd_edit', 're_edit', 'rf_edit'};

for i = 1:length(control_fields)
    if isfield(handles, control_fields{i})
        set(handles.(control_fields{i}), 'Enable', 'on');
        set(handles.(control_fields{i}), 'BackgroundColor', 'white');
    end
end
end

function enable_calculation_controls(handles)
% Enable calculation controls for non-SCN materials
calc_control_fields = {'dt_min_edit', 'dt_max_edit', 'sampling_interval_edit', 'split_pt_edit', ...
                       'lower_density_mult_edit', 'upper_density_mult_edit', 'low_density_edit', 'high_density_edit', ...
                       'c0_values_edit', 'v_a_edit', 'v_b_edit', 'v_c_edit', 'rd_edit', 're_edit', 'rf_edit'};

for i = 1:length(calc_control_fields)
    if isfield(handles, calc_control_fields{i})
        set(handles.(calc_control_fields{i}), 'Enable', 'on');
        set(handles.(calc_control_fields{i}), 'BackgroundColor', 'white');
    end
end
end

function disable_material_parameters(handles)
% Disable material parameter controls for predefined materials
material_control_fields = {'cp_edit', 'dh_edit', 'a_edit', 'gibbs_edit', 'dal_edit', 'mal_edit', 'k0al_edit', 'sigma_edit'};

for i = 1:length(material_control_fields)
    if isfield(handles, material_control_fields{i})
        set(handles.(material_control_fields{i}), 'Enable', 'off');
        set(handles.(material_control_fields{i}), 'BackgroundColor', [0.9, 0.9, 0.9]);
    end
end
end


function test_parameter_save_load()
% Test function to verify parameter save/load functionality
fprintf('=== TESTING PARAMETER SAVE/LOAD ===\n');

% Check if file exists
if exist('material_params.mat', 'file')
    fprintf('File exists: material_params.mat\n');
    
    % Load and display contents
    data = load('material_params.mat');
    if isfield(data, 'params')
        fprintf('Params structure found in file\n');
        fprintf('Available fields:\n');
        fields = fieldnames(data.params);
        for i = 1:length(fields)
            field_name = fields{i};
            field_value = data.params.(field_name);
            if isnumeric(field_value)
                fprintf('  %s: %g\n', field_name, field_value);
            else
                fprintf('  %s: %s\n', field_name, field_value);
            end
        end
    else
        fprintf('ERROR: No params structure found in file\n');
    end
else
    fprintf('File does NOT exist: material_params.mat\n');
end

fprintf('Current working directory: %s\n', pwd);
fprintf('Files in current directory:\n');
files = dir('*.mat');
for i = 1:length(files)
    fprintf('  %s\n', files(i).name);
end
fprintf('================================\n');
end


function cleanup_parallel_pool_on_close()
% CLEANUP_PARALLEL_POOL_ON_CLOSE - Clean up parallel pool when GUI closes
try
    pool = gcp('nocreate');
    if ~isempty(pool)
        fprintf('Closing parallel pool (%d workers) as GUI is closing...\n', pool.NumWorkers);
        delete(pool);
        fprintf('✓ Parallel pool closed successfully\n');
    else
        fprintf('No parallel pool to close\n');
    end
catch ME
    fprintf('⚠ Warning: Error closing parallel pool: %s\n', ME.message);
end
end

function initialize_parallel_pool_for_calculation()
% INITIALIZE_PARALLEL_POOL_FOR_CALCULATION - Initialize parallel pool when starting calculation
% This function ensures parallel pool is ready before calculation starts

fprintf('=== Initializing Parallel Computing Pool for Calculation ===\n');

% Check if parallel pool already exists
pool = gcp('nocreate');
if ~isempty(pool)
    fprintf('✓ Parallel pool already exists (%d workers) - ready for calculation\n', pool.NumWorkers);
    return;
end

% Check if Parallel Computing Toolbox is available
if ~license('test', 'Distrib_Computing_Toolbox')
    fprintf('⚠ Parallel Computing Toolbox not available - will use serial computation\n');
    return;
end

% Determine optimal number of workers
num_physical_cores = feature('numcores');
num_logical_cores = java.lang.Runtime.getRuntime.availableProcessors;
optimal_workers = max(2, min(num_physical_cores, num_logical_cores - 1));

fprintf('Creating parallel pool with %d workers for calculation...\n', optimal_workers);

try
    % Create the pool synchronously (blocking)
    pool = parpool('local', optimal_workers);
    
    fprintf('✓ Parallel pool created successfully (%d workers)\n', pool.NumWorkers);
    fprintf('  Pool will remain active until GUI is closed\n');
    
catch ME
    fprintf('⚠ Failed to create parallel pool: %s\n', ME.message);
    fprintf('  Will use serial computation for this calculation\n');
end
end

