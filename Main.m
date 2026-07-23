% This function calls all code and generates all figures from the Ron et al
% manuscript.

% TODO - add instructions for how to install SSIT

% Add source code to path.
addpath(genpath('src'));

opts.GRModels = {'all'}; 
% Select 'all', 'none' or subset of {'DexControlsGRImport_NucCytDeg',
%   'DexControlsGRImport_NucDeg','DexControlsGRImport_CytDeg',
%   'DexControlsGRExport_NucCytDeg'}

opts.FullModels = {'all'};
% Select 'all', 'none', or subset of {'TVandSat',
%   'NoTV','NoSat','LinDecay','Mechanistic'}

opts.regenLibrary = false;% Regenerate model library.  Usually only needed in subsequent runs when switching between machines. 
opts.remakeFigs = false;   % Remake figures (true) or use saved (false)? NOTE: Can take half an hour or more.
opts.refitModels = false; % Refit models (true) or use saved (false)?   NOTE: This can take many minutes to hours per model.
opts.hideFigures = true;  % Hide annoying figures during generation.
opts.recomputeSSA = false; % Re-run a new set of SSA for each prediction of figure to be generated (full models only).

% Call function to generate all final plots
Wrapper_For_Final_Plots(opts);

disp('********************************************************')
disp('Writing Table of Parameter Comparisons')
disp('********************************************************')
DisplayParameters


disp('********************************************************')
disp('Generating Final Nuc/Cyt Comparison Plots')
disp('********************************************************')

% Call function to generate the plots of means and variaces versus time.
Wrapper_For_Final_Model_Comparisons(opts);

disp('********************************************************')
disp('Generating Verification of Hybrid FSP against SSA')
disp('********************************************************')
% Call function to solve Hybrid FSP model and verify against SSA.
verifyHybridFSP
