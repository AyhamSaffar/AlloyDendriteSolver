%% LGK_model.m - Lipton-Glicksman-Kurz Dendritic Growth ModelAuthor: Xinyi Hao
% Author: Xinyi Hao
% Based on core code from Sihan Sun
% ========================================================================
% 
% PROGRAMME OVERVIEW:
% This MATLAB programme suite implements the Lipton-Glicksman-Kurz (LGK) 
% model for free dendritic growth in undercooled binary alloy melts. The 
% model calculates dendrite tip velocity (V) and tip radius (R) by coupling 
% heat and mass transport around a parabolic dendrite tip with morphological 
% stability criteria.
%
% PHYSICAL BASIS:
% Free dendrites growing in liquid alloys reject both latent heat and 
% solute (when k₀ ≠ 1). The LGK model combines:
% 1. Ivantsov's transport solution for diffusion fields around parabolic tips
% 2. Marginal stability criterion determining the operating tip radius
% 3. Coupled thermal and solutal undercooling balance equations
%
% The model predicts that:
% - Dendrite tip radius passes through a minimum with increasing solute concentration
% - Growth velocity increases with solute concentration then decreases at higher levels
% - Both phenomena result from the competition between thermal and solutal effects
%
% MATHEMATICAL FORMULATION:
% The model solves two coupled nonlinear equations:
% 1. Undercooling balance: ΔT = ΔTₜ + ΔTc + ΔTᵣ
%    where ΔTₜ = thermal, ΔTc = solutal, ΔTᵣ = curvature undercooling
% 2. Stability criterion: R = σ*Γ/(mGc - G)
%    where σ is stability constant, Γ is Gibbs-Thomson coefficient
%
% PROGRAMME STRUCTURE:
% - MaterialControlPanel.m: Main GUI for parameter input and calculation control
% - calculate_VR.m: Core Newton-Raphson solver for V and R values
% - run_main_calculation.m: Orchestrates multi-concentration calculations
% - create_deltaT_sequence.m: Generates undercooling value sequences
% - create_analytical_guess_functions.m: Provides initial guess functions
%
% MATERIALS SUPPORTED:
% Pre-configured parameters for: AZ91, Al-4wt%Cu, Al-Cu, Al-Fe, Sn-Ag, 
% Sn-Cu, SCN-Acetone, Mg-Alloy, plus custom material definition capability
%
% CALCULATION FEATURES:
% - Parallel processing support for multiple concentrations
% - Real-time progress monitoring with GUI feedback
% - Automatic parameter fitting and optimisation
% - Comprehensive result visualisation and Excel export
% - Both parametric and analytical initial guess functions
%
% REFERENCE:
% Lipton, J., Glicksman, M.E., Kurz, W. (1984). "Dendritic Growth into 
% Undercooled Alloy Melts". Materials Science and Engineering, 65, 57-63.
%
% AUTHORS: Based on LGK theory, implemented with modern MATLAB features
% VERSION: Advanced GUI version with parallel processing capabilities
% ========================================================================

% Clear workspace and close existing figures to ensure clean start
clear all
close all
clc

% Launch the material control panel GUI - this is the main interface
% The GUI handles all user interactions, parameter input, and calculation management
fprintf('Starting Material Parameter Control Panel...\n');
MaterialControlPanel;

% The GUI will handle the rest of the calculation workflow
% The original main calculation code is now managed through the GUI interface
% This prevents conflicts when the GUI is active and ensures proper resource management
return;