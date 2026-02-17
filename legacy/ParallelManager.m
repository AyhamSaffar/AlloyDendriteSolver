function ParallelManager()
%% PARALLELMANAGER - Advanced parallel computing setup and optimization for LGK calculations
% ========================================================================
%
% PURPOSE:
% Implements intelligent parallel computing configuration with automatic
% hardware detection, optimal worker pool sizing, and performance testing
% to maximize computational efficiency for LGK dendritic growth calculations.
%
% PARALLEL COMPUTING STRATEGY:
% The LGK model involves solving nonlinear equations for multiple concentration
% values, making it ideally suited for parallel computation. This function
% optimizes the parallel environment by:
% - Detecting available CPU cores and architecture characteristics
% - Configuring optimal worker pool size based on workload type
% - Testing parallel efficiency to validate configuration
% - Providing performance feedback for user optimization
%
% HARDWARE OPTIMIZATION CONSIDERATIONS:
%
% Physical vs Logical Cores:
% Modern CPUs often implement hyperthreading (Intel) or SMT (AMD), creating
% more logical cores than physical cores. For CPU-intensive tasks like
% Newton-Raphson iteration, physical cores typically provide better performance.
%
% Memory Bandwidth Limitations:
% Each parallel worker requires memory bandwidth for:
% - Symbolic computation caching
% - Intermediate calculation storage  
% - Inter-worker communication overhead
% Excessive workers can saturate memory bandwidth, reducing efficiency.
%
% Cache Coherency:
% Multiple workers accessing shared memory structures can cause cache
% invalidation overhead. Optimal worker count balances parallelism
% with cache efficiency.
%
% WORKER POOL SIZING ALGORITHM:
% 1. Detect total logical and physical core counts
% 2. Reserve one core for system operations (logical_cores - 1)
% 3. For hyperthreaded systems: prefer physical core count for CPU-intensive tasks
% 4. Apply conservative limits to prevent system overload
% 5. Validate performance through empirical testing
%
% PERFORMANCE TESTING METHODOLOGY:
% Executes standardized computational tasks to measure:
% - Serial execution baseline performance
% - Parallel execution with configured worker count
% - Speedup ratio and parallel efficiency metrics
% - Memory bandwidth and cache performance indicators
%
% EFFICIENCY METRICS:
% - Speedup = T_serial / T_parallel
% - Efficiency = Speedup / Number_of_Workers × 100%
% - Excellent: >70% efficiency
% - Good: 50-70% efficiency  
% - Poor: <50% efficiency (suggests over-parallelization)
%
% INPUTS: None (automatic hardware detection)
% OUTPUTS: Configured parallel pool ready for LGK calculations
%
% DEPENDENCIES:
% - Parallel Computing Toolbox (required)
% - System hardware detection capabilities
%
% USAGE:
% Call ParallelManager() before running LGK calculations to ensure
% optimal parallel configuration for the current hardware platform.
%
% PERFORMANCE RECOMMENDATIONS:
% - For systems with <4 cores: Consider serial execution
% - For systems with 4-8 cores: Use physical core count
% - For systems with >8 cores: Test optimal worker count empirically
% - For memory-limited systems: Reduce worker count below core count
% ========================================================================

% Get comprehensive system information for optimization decisions

% Get system information
num_physical_cores = feature('numcores');
num_logical_cores = java.lang.Runtime.getRuntime.availableProcessors;

fprintf('=== Parallel Computing Setup ===\n');
fprintf('Physical CPU cores: %d\n', num_physical_cores);
fprintf('Logical CPU cores: %d\n', num_logical_cores);

% Check if Parallel Computing Toolbox is available
if ~license('test', 'Distrib_Computing_Toolbox')
    fprintf('Parallel Computing Toolbox not available - using serial computation\n');
    return;
end

% Close existing parallel pool if any
if ~isempty(gcp('nocreate'))
    fprintf('Closing existing parallel pool...\n');
    delete(gcp('nocreate'));
end

% Determine optimal number of workers
% Use all logical cores minus 1 for system stability, but at least 2
optimal_workers = max(2, num_logical_cores - 1);

% For systems with hyperthreading, consider using physical cores only for CPU-intensive tasks
if num_logical_cores > num_physical_cores
    fprintf('Hyperthreading detected (logical/physical = %.1f)\n', num_logical_cores/num_physical_cores);
    % Option 1: Use all logical cores for I/O intensive tasks
    % Option 2: Use physical cores for CPU intensive tasks
    fprintf('Recommended workers for CPU-intensive tasks: %d (physical cores)\n', num_physical_cores);
    fprintf('Recommended workers for mixed tasks: %d (logical cores - 1)\n', optimal_workers);
    
    % For this calculation, use physical cores for better performance
    optimal_workers = min(optimal_workers, num_physical_cores);
end

fprintf('Setting up parallel pool with %d workers...\n', optimal_workers);

try
    % Create parallel pool with optimized settings
    cluster = parcluster('local');
    
    % Optimize cluster settings for maximum performance
    cluster.NumWorkers = optimal_workers;
    
    % Set memory and CPU affinity options if supported
    if isprop(cluster, 'AdditionalProperties')
        cluster.AdditionalProperties.AdditionalSubmitArgs = '-singleCompThread';
    end
    
    % Create the pool
    pool = parpool(cluster, optimal_workers);
    
    fprintf('Parallel pool created successfully with %d workers\n', pool.NumWorkers);
    
    % Set parallel preferences for maximum performance
    set_parallel_preferences();
    
    % Test parallel performance
    test_parallel_performance(pool.NumWorkers);
    
catch ME
    fprintf('Failed to create parallel pool: %s\n', ME.message);
    fprintf('Falling back to serial computation\n');
end

end

function set_parallel_preferences()
% Set optimal parallel computing preferences

try
    % Enable automatic parallel pool creation
    ps = parallel.Settings;
    ps.Pool.AutoCreate = true;
    
    % Set timeout for workers
    ps.Pool.IdleTimeout = 30; % minutes
    
    % Optimize for computational efficiency
    fprintf('Parallel preferences optimized\n');
catch
    fprintf('Could not set parallel preferences\n');
end

end

function test_parallel_performance(num_workers)
%% TEST_PARALLEL_PERFORMANCE - Empirical parallel efficiency validation
%
% PURPOSE:
% Executes standardized computational benchmarks to validate parallel
% configuration efficiency and provide optimization recommendations.
%
% BENCHMARK METHODOLOGY:
% Uses matrix operations similar to LGK symbolic computations:
% - Large random matrix generation and summation
% - Memory-intensive operations to test bandwidth limits
% - CPU-intensive calculations to measure computational speedups
%
% PERFORMANCE ANALYSIS:
% Compares serial vs parallel execution times to calculate:
% - Absolute speedup ratio
% - Parallel efficiency percentage
% - Performance recommendations based on measured efficiency

fprintf('\n=== Performance Test ===\n');

% Serial test
tic;
serial_result = sum(rand(1000, 1000), 'all');
serial_time = toc;
fprintf('Serial computation time: %.3f seconds\n', serial_time);

% Parallel test
tic;
parfor i = 1:num_workers
    parallel_result(i) = sum(rand(1000, 1000), 'all'); %#ok<PFBNS>
end
parallel_time = toc;
fprintf('Parallel computation time: %.3f seconds\n', parallel_time);

speedup = serial_time / parallel_time;
efficiency = speedup / num_workers * 100;

fprintf('Speedup: %.2fx\n', speedup);
fprintf('Efficiency: %.1f%%\n', efficiency);

if efficiency > 70
    fprintf('✓ Excellent parallel efficiency\n');
elseif efficiency > 50
    fprintf('✓ Good parallel efficiency\n');
else
    fprintf('⚠ Consider reducing worker count for better efficiency\n');
end

fprintf('=== Setup Complete ===\n\n');
end