function deltaT_values = create_deltaT_sequence(varargin)
%% CREATE_DELTAT_SEQUENCE - Generate adaptive undercooling sequences for LGK calculations
% ========================================================================
%
% PURPOSE:
% Creates sophisticated undercooling (ΔT) sequences with adaptive density
% control for optimal computational efficiency in LGK dendritic growth
% calculations. The function supports both single-region linear sequences
% and two-region interpolation with density transformations.
%
% SCIENTIFIC RATIONALE:
% Dendritic growth exhibits different physical regimes across the undercooling
% range, requiring non-uniform sampling for computational efficiency:
% - Low ΔT: Dominated by surface tension effects, requires fine resolution
% - High ΔT: Dominated by kinetic effects, can use coarser sampling
% - Transition regions: May require enhanced density for regime changes
%
% TWO-MODE OPERATION:
%
% MODE 1: Simple Linear Sequence
% Used when split_point is invalid (NaN, outside range, or disabled)
% Creates uniform spacing: deltaT_min : sampling_interval : deltaT_max
% Density control parameters are ignored in this mode
%
% MODE 2: Two-Region Interpolation with Density Control
% Used when split_point is valid (within deltaT_min < split_point < deltaT_max)
% Creates two regions with independent sampling characteristics:
% - Lower region: [deltaT_min, split_point] with enhanced density control
% - Upper region: [split_point, deltaT_max] with independent density control
%
% DENSITY TRANSFORMATION ALGORITHM:
% Applies power-law transformations to achieve non-uniform point distribution:
%
% For density > 1: t_transformed = t^(1/density)
% - Concentrates points toward the start of the region
% - Useful for capturing rapid changes near critical points
%
% For density < 1: t_transformed = 1 - (1-t)^(1/density)  
% - Concentates points toward the end of the region
% - Useful for emphasizing high-undercooling behavior
%
% For density = 1: Linear spacing (no transformation)
%
% INPUTS (Name-Value Pairs):
%   'deltaT_min' - Minimum undercooling value [K] (default: 0.1)
%   'deltaT_max' - Maximum undercooling value [K] (default: 20)
%   'sampling_interval' - Base sampling interval [K] (default: 0.5)
%   'split_point' - Transition point between regions [K] (default: NaN)
%   'lower_density_mult' - Lower region density multiplier [-] (default: 5.0)
%   'upper_density_mult' - Upper region density multiplier [-] (default: 1.0)  
%   'low_end_density' - Low end density transformation factor [-] (default: 3.0)
%   'high_end_density' - High end density transformation factor [-] (default: 1.0)
%
% DENSITY MULTIPLIER INTERPRETATION:
% - lower_density_mult = 5.0: 5× more points in lower region than base interval
% - upper_density_mult = 1.0: Standard density in upper region
% - Combined effect: Lower region gets 5× finer sampling than upper region
%
% DENSITY TRANSFORMATION INTERPRETATION:
% - low_end_density > 1: More points near deltaT_min
% - low_end_density < 1: More points near split_point
% - high_end_density > 1: More points near deltaT_max  
% - high_end_density < 1: More points near split_point
%
% OUTPUTS:
%   deltaT_values - Vector of undercooling values in descending order [K]
%                   Optimized for Newton-Raphson convergence (high to low ΔT)
%
% ALGORITHM WORKFLOW:
% 1. Parse and validate input parameters
% 2. Determine operation mode based on split_point validity
% 3. Calculate region-specific intervals and point counts
% 4. Generate base sequences for each region
% 5. Apply density transformations to interior points
% 6. Combine regions and remove duplicates
% 7. Sort in descending order for optimal solver performance
% 8. Provide comprehensive diagnostic output
%
% USAGE EXAMPLES:
%
% Simple linear sequence:
%   seq = create_deltaT_sequence('deltaT_min', 0.5, 'deltaT_max', 10, 'sampling_interval', 0.2);
%
% Two-region with enhanced lower density:
%   seq = create_deltaT_sequence('deltaT_min', 0.1, 'deltaT_max', 20, ...
%                                'split_point', 2.0, 'lower_density_mult', 10, ...
%                                'low_end_density', 2.0);
%
% COMPUTATIONAL ADVANTAGES:
% - Adaptive sampling reduces total calculation points
% - Enhanced resolution in physically critical regions
% - Smooth transitions between different sampling densities
% - Eliminates redundant calculations in over-sampled regions
%
% PERFORMANCE CONSIDERATIONS:
% - Descending order optimizes Newton-Raphson initial guess inheritance
% - Duplicate removal prevents redundant calculations
% - Diagnostic output aids in sequence optimization
% - Memory-efficient generation for large sequences
% ========================================================================

% Parse input arguments
p = inputParser;
addParameter(p, 'deltaT_min', 0.1, @(x) isnumeric(x) && x > 0);
addParameter(p, 'deltaT_max', 20, @(x) isnumeric(x) && x > 0);
addParameter(p, 'sampling_interval', 0.5, @(x) isnumeric(x) && x > 0);
addParameter(p, 'split_point', NaN, @(x) isnumeric(x));
addParameter(p, 'lower_density_mult', 5.0, @(x) isnumeric(x) && x > 0);
addParameter(p, 'upper_density_mult', 1.0, @(x) isnumeric(x) && x > 0);
addParameter(p, 'low_end_density', 3.0, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'high_end_density', 1.0, @(x) isnumeric(x) && x >= 0);
parse(p, varargin{:});

deltaT_min = p.Results.deltaT_min;
deltaT_max = p.Results.deltaT_max;
sampling_interval = p.Results.sampling_interval;
split_point = p.Results.split_point;
lower_density_mult = p.Results.lower_density_mult;
upper_density_mult = p.Results.upper_density_mult;
low_end_density = p.Results.low_end_density;
high_end_density = p.Results.high_end_density;

fprintf('=== deltaT Sequence Generation ===\n');
fprintf('Input parameters:\n');
fprintf('  deltaT_min: %.3f, deltaT_max: %.3f\n', deltaT_min, deltaT_max);
fprintf('  sampling_interval: %.3f\n', sampling_interval);
fprintf('  split_point: %.3f\n', split_point);
fprintf('  lower_density_mult: %.2f, upper_density_mult: %.2f\n', lower_density_mult, upper_density_mult);
fprintf('  low_end_density: %.2f, high_end_density: %.2f\n', low_end_density, high_end_density);

% Validate input
if deltaT_min >= deltaT_max
    error('deltaT_min must be less than deltaT_max');
end

% Check if split point is valid (within range, exclusive)
is_split_valid = ~isnan(split_point) && split_point > deltaT_min && split_point < deltaT_max;

if ~is_split_valid
    % Mode 1: Simple linear sequence (density controls ignored)
    fprintf('\nMode 1: Split point invalid or outside range\n');
    fprintf('Using simple sequence: deltaT_min:sampling_interval:deltaT_max\n');
    
    % Create simple linear sequence
    deltaT_values = deltaT_min:sampling_interval:deltaT_max;
    
    % Ensure deltaT_max is included if not exactly reached
    if abs(deltaT_values(end) - deltaT_max) > 1e-10
        deltaT_values = [deltaT_values, deltaT_max];
    end
    
    % Sort in descending order
    deltaT_values = sort(deltaT_values, 'descend');
    
    fprintf('Generated %d points from %.3f to %.3f\n', length(deltaT_values), min(deltaT_values), max(deltaT_values));
    
else
    % Mode 2: Two-region interpolation with density control
    fprintf('\nMode 2: Two-region interpolation (split point valid)\n');
    fprintf('Split point %.3f is valid (%.3f < %.3f < %.3f)\n', split_point, deltaT_min, split_point, deltaT_max);
    
    % Calculate intervals for each region
    lower_interval = sampling_interval / lower_density_mult;
    upper_interval = sampling_interval / upper_density_mult;
    
    fprintf('Region intervals:\n');
    fprintf('  Lower region [%.3f, %.3f]: interval = %.6f\n', deltaT_min, split_point, lower_interval);
    fprintf('  Upper region [%.3f, %.3f]: interval = %.6f\n', split_point, deltaT_max, upper_interval);
    
    fprintf('Applying density transformations:\n');
    fprintf('  Low end density: %.2f (>1: more points near deltaT_min, <1: more points near split_point)\n', low_end_density);
    fprintf('  High end density: %.2f (>1: more points near deltaT_max, <1: more points near split_point)\n', high_end_density);

    % Generate lower region sequence [deltaT_min, split_point]
    lower_seq = deltaT_min:lower_interval:split_point;
    if abs(lower_seq(end) - split_point) > 1e-10
        lower_seq = [lower_seq, split_point];
    end
    
    % Apply density transformation to lower region (excluding endpoints)
    if length(lower_seq) > 2 && low_end_density ~= 1.0
        % Transform interior points
        interior_points = lower_seq(2:end-1);
        if ~isempty(interior_points)
            % Normalize to [0,1] where 0=deltaT_min, 1=split_point
            t_norm = (interior_points - deltaT_min) / (split_point - deltaT_min);
            
            % Apply density transformation:
            % low_end_density > 1: more points near deltaT_min (t=0)
            % low_end_density < 1: more points near split_point (t=1)
            if low_end_density > 1
                % More density near deltaT_min: use t^(1/density)
                t_transformed = t_norm.^(1/low_end_density);
            else
                % More density near split_point: use 1-(1-t)^(1/density)
                t_transformed = 1 - (1 - t_norm).^(1/low_end_density);
            end
            
            interior_transformed = deltaT_min + t_transformed * (split_point - deltaT_min);
            lower_seq = [deltaT_min, interior_transformed, split_point];
        end
    end
    
    % Generate upper region sequence [split_point, deltaT_max]
    upper_seq = split_point:upper_interval:deltaT_max;
    if abs(upper_seq(end) - deltaT_max) > 1e-10
        upper_seq = [upper_seq, deltaT_max];
    end
    
    % Apply density transformation to upper region (excluding endpoints)
    if length(upper_seq) > 2 && high_end_density ~= 1.0
        % Transform interior points
        interior_points = upper_seq(2:end-1);
        if ~isempty(interior_points)
            % Normalize to [0,1] where 0=split_point, 1=deltaT_max
            t_norm = (interior_points - split_point) / (deltaT_max - split_point);
            
            % Apply density transformation:
            % high_end_density > 1: more points near deltaT_max (t=1)
            % high_end_density < 1: more points near split_point (t=0)
            if high_end_density > 1
                % More density near deltaT_max: use 1-(1-t)^(1/density)
                t_transformed = 1 - (1 - t_norm).^(1/high_end_density);
            else
                % More density near split_point: use t^(1/density)
                t_transformed = t_norm.^(1/high_end_density);
            end
            
            interior_transformed = split_point + t_transformed * (deltaT_max - split_point);
            upper_seq = [split_point, interior_transformed, deltaT_max];
        end
    end
    
    % Remove split_point from one sequence to avoid duplication
    upper_seq = upper_seq(2:end);  % Remove first element (split_point)
    
    % Combine sequences
    deltaT_values = [lower_seq, upper_seq];
    
    % Sort in descending order
    deltaT_values = sort(deltaT_values, 'descend');
    
    % Remove duplicates with small tolerance
    deltaT_values = remove_duplicates(deltaT_values, 1e-10);
    
    fprintf('Region results:\n');
    fprintf('  Lower region: %d points\n', length(lower_seq));
    fprintf('  Upper region: %d points\n', length(upper_seq));
    fprintf('  Combined: %d points\n', length(deltaT_values));
end

% Final validation and output
fprintf('\nFinal sequence:\n');
fprintf('  Total points: %d\n', length(deltaT_values));
fprintf('  Range: %.6f to %.6f\n', min(deltaT_values), max(deltaT_values));

% Show first and last few points
fprintf('  First 5 points: ');
for i = 1:min(5, length(deltaT_values))
    fprintf('%.4f ', deltaT_values(i));
end
fprintf('\n  Last 5 points: ');
for i = max(1, length(deltaT_values)-4):length(deltaT_values)
    fprintf('%.4f ', deltaT_values(i));
end
fprintf('\n');

% Check spacing statistics
if length(deltaT_values) > 1
    spacings = abs(diff(deltaT_values));
    fprintf('  Spacing - Min: %.6f, Max: %.6f, Avg: %.6f\n', ...
            min(spacings), max(spacings), mean(spacings));
end

fprintf('=== End deltaT Sequence Generation ===\n\n');

end

function clean_values = remove_duplicates(values, tolerance)
% Remove values that are too close to each other
if length(values) <= 1
    clean_values = values;
    return;
end

clean_values = values(1);
for i = 2:length(values)
    if abs(values(i) - clean_values(end)) > tolerance
        clean_values = [clean_values, values(i)];
    end
end
end
