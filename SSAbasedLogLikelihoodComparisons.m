%% Testing Approximate Log-Likelihood calculations for SSA.
% Run from the repository root. Requires SSIT on the path, e.g.,
%   addpath(genpath('src')); addpath(genpath('../../SSIT'));
addpath(genpath('src'));
addpath(genpath('tmpPropensityFunctions'));

% Choose what to run.an
runCalibration = false;    % SSA likelihood noise vs number of runs (results
                          % already in savedResults/SSALikelihood).
runSeqABC = true;         % Sequential ABC (MH) for the TVandSat model.

%% SSA likelihood calibration for the full TVandSat model (DUSP1 model 2)
% Mean and std of the SSA/KDE log-likelihood of the joint nuclear and
% cytoplasmic DUSP1 mRNA data (100 nM Dex) versus the number of SSA runs,
% at the final fitted parameters. The full model is too large for the FSP,
% so there is no exact reference value here.
if runCalibration
    calTVandSat = ssaCalibrationTVandSat(nSimsList=[2500, 10000, 40000], nReps=8);
end

%% Sequential ABC for the full TVandSat model
% Metropolis-Hastings with the SSA/KDE likelihood of the nuclear and
% cytoplasmic DUSP1 data (SSIT.runSequentialABC), plus the (FSP-based)
% transcription-site likelihood and the log-normal prior, with the MH
% settings and proposal tuning of GR_DUSP1_Fitting_Final (model 2).
% Parameters are read from and saved to
% savedParameters/Parameters_Final_DUSP1_w_GR_Binning (created from
% Parameters_Final, which is never written); chains and calibration go
% to savedResults/SSALikelihood.
if runSeqABC
    abcTVandSat = seqABCTVandSat(marginalsOnly=false);
end

%% Local functions
function results = ssaCalibrationTVandSat(opts)
% ssaCalibrationTVandSat - Calibrate the SSA-based approximate
% log-likelihood for the full TVandSat DUSP1 model (time-varying GR input
% and saturating degradation; model 2 in GR_DUSP1_Fitting_Final).
%
%   results = ssaCalibrationTVandSat()
%   results = ssaCalibrationTVandSat(nSimsList=[100 1000], nReps=10, ...)
%
% Steps:
%   1. Load 'fullSSAModel_100_2StateGen' from the model library. Its data
%      are the joint nuclear (rna <- num_nuc_noTS) and cytoplasmic
%      (rCyt <- num_cyto_spots) DUSP1 counts at 100 nM Dex.
%   2. Set the final parameters exactly as GR_DUSP1_Fitting_Final does for
%      model 2: GR parameters (GRpars), nuclear DUSP1 parameters
%      (Pars_KON_KOFF_KR_zeroOFF), parameter 23 set to zero, and the
%      fitted free parameters (FinalModel_allButGR).
%   3. Run SSIT.ssaLikelihoodConvergence, which evaluates the SSA-only
%      objective (default 'kde_loglik' on the JOINT nuclear/cytoplasmic
%      data) nReps times for each number of SSA runs per time point.
%
% Options:
%   nSimsList         -- SSA runs per time point ([10 30 100 300 1000]).
%   nReps             -- repetitions per number of runs (5).
%   lossFun           -- objective for computeLossFunctionSSA ('kde_loglik').
%   kdeBiasCorrection -- 'jensen' (default) or 'none'.
%   kdeBandwidth      -- KDE bandwidth rule ('lcv').
%   parsFile          -- parameter file ('savedParameters/Parameters_Final').
%   modelLibrary      -- model library ('savedParameters/GRDusp1ModelLibrary').
%   saveResults       -- (true) save results and figure under
%                        savedResults/SSALikelihood and
%                        savedFigures/SSALikelihood (new folders; nothing
%                        existing is overwritten except earlier output of
%                        this function with the same tag).
%   tag               -- ('TVandSat') name used in the saved files.
%
% Output: the struct from ssaLikelihoodConvergence (nSims, values, mean,
% std, secondsPerEval, figureHandle, ...) plus the model, the parameter
% file and the options used.
arguments
    opts.nSimsList (1,:) double {mustBePositive,mustBeInteger} = [2500, 10000, 40000]
    opts.nReps (1,1) double {mustBePositive,mustBeInteger} = 5
    opts.lossFun = 'kde_loglik'
    opts.kdeBiasCorrection {mustBeMember(opts.kdeBiasCorrection,{'none','jensen'})} = 'jensen'
    opts.kdeBandwidth = 'lcv'
    opts.parsFile = 'savedParameters/Parameters_Final'
    opts.modelLibrary = 'savedParameters/GRDusp1ModelLibrary'
    opts.saveResults (1,1) logical = true
    opts.tag = 'TVandSat'
end

% 1. Model with joint nuclear/cytoplasmic DUSP1 data at 100 nM Dex.
SSAModel = dusp1ModelLibrary_Final('fullSSAModel_100_2StateGen', false, opts.modelLibrary);

% 2. Final parameters (as in GR_DUSP1_Fitting_Final, model 2).
indsFreePars = [1,2,3,4,14,15,16,19,20,21,22,24,25,26];
zeroPars = 23;
P = load(opts.parsFile, 'GRpars', 'Pars_KON_KOFF_KR_zeroOFF', 'FinalModel_allButGR');
SSAModel.parameters(5:12,2) = num2cell(P.GRpars);
SSAModel.parameters([1,2,3,4,19,24,25,26],2) = num2cell(P.Pars_KON_KOFF_KR_zeroOFF);
SSAModel.parameters(zeroPars,2) = num2cell(zeros(size(zeroPars)));
SSAModel.parameters(indsFreePars,2) = num2cell(P.FinalModel_allButGR);
SSAModel.fittingOptions.modelVarsToFit = indsFreePars;
SSAModel.solutionScheme = 'SSA';
if isunix && ~ismac
    % Same workaround as GR_DUSP1_Fitting_Final for nonlinear propensities.
    SSAModel.ssaOptions.useC = false;
end
SSAModel=SSAModel.solve(solver='ssa');
fprintf('TVandSat model: %s fitted to species {%s}; %d data times; %d cells.\n', ...
    opts.parsFile, strjoin(SSAModel.dataSet.linkedSpecies.keys(), ', '), ...
    numel(SSAModel.dataSet.times), sum(SSAModel.dataSet.nCells));

% 3. SSA-only likelihood calibration (no FSP reference: the full model is
% too large for the FSP).
results = SSAModel.ssaLikelihoodConvergence(nSimsList=opts.nSimsList, ...
    nReps=opts.nReps, lossFun=opts.lossFun, compareFSP=false, makePlot=true, ...
    kdeBandwidth=opts.kdeBandwidth, kdeBiasCorrection=opts.kdeBiasCorrection);
results.model = SSAModel;
results.parsFile = opts.parsFile;
results.options = opts;
if isfield(results, 'figureHandle')
    title(results.figureHandle.Children(1), sprintf(['TVandSat model, joint nuclear/cytoplasmic DUSP1 (100 nM Dex): ', ...
        'SSA objective ''%s'' (correction: %s), %d repeats per point'], ...
        opts.lossFun, opts.kdeBiasCorrection, opts.nReps), 'Interpreter', 'none');
end

if opts.saveResults
    resDir = 'savedResults/SSALikelihood';
    figDir = 'savedFigures/SSALikelihood';
    if ~isfolder(resDir), mkdir(resDir); end
    if ~isfolder(figDir), mkdir(figDir); end
    baseName = sprintf('SSALikelihoodCalibration_%s_%s', opts.tag, opts.kdeBiasCorrection);
    toSave = rmfield(results, intersect(fieldnames(results), {'figureHandle','model'}));
    save(fullfile(resDir, baseName), '-struct', 'toSave');
    if isfield(results, 'figureHandle')
        savefig(results.figureHandle, fullfile(figDir, [baseName, '.fig']));
        exportgraphics(results.figureHandle, fullfile(figDir, [baseName, '.png']), 'Resolution', 130);
    end
    fprintf('Saved %s and figure in %s.\n', fullfile(resDir, baseName), figDir);
end
end

function out = seqABCTVandSat(opts)
% seqABCTVandSat - Sequential ABC (SSIT.runSequentialABC) for the full
% TVandSat DUSP1 model (model 2 of GR_DUSP1_Fitting_Final), with the MH
% settings and proposal tuning of GR_DUSP1_Fitting_Final.
%
%   out = seqABCTVandSat()
%   out = seqABCTVandSat(marginalsOnly=true, finalSamples=2000, ...)
%
% Target (log posterior, natural log):
%   log L_SSA(nuclear/cytoplasmic DUSP1 data)   -- SSA/KDE estimate, with
%                                                 the sequential test
% + log L_TS(transcription-site data)          -- FSP (computeTSlikelihood)
% + log prior                                   -- as in computeTotalObjFun:
%       -sum((log10(x)-log10PriorMean).^2 ./ (2*log10PriorStd.^2))
% The TS likelihood and prior are deterministic, so they enter the MH
% threshold (runSequentialABC's logPrior input) and only the SSA term is
% estimated sequentially. Unlike GR_DUSP1_Fitting_Final, the three terms
% are not weighted (weights_ABC = [250 47.8 1] there): all three are now
% log-likelihoods/log-priors, so they are added directly.
%
% Proposal tuning (as in GR_DUSP1_Fitting_Final):
%   * proposal x + N(0, scaleMH*covPost + 1e-6*I) in log-parameter space,
%     starting from covPost_<parName> in parsFile (or 1e-4*I if absent)
%     and scaleMH = 1;
%   * tuning runs of tuneSamples samples (thin tuneThin); after each run,
%     scaleMH is multiplied by 0.2 (acceptance < 0.05), 0.75 (< 0.2),
%     3 (> 0.6) or 1.5 (> 0.45); otherwise the final run (finalSamples,
%     thin finalThin) starts;
%   * after each run, if the number of accepted samples exceeds the
%     number of parameters, covPost = cov(samples) and the starting point
%     is the best parameters found.
%
% SSA runs: stages of levels = [2500 5000 10000] runs per time point
% (with enforceIndependence = false, a 2500-run evaluation costs about
% the same as a 300-run one because of fixed overheads). With
% enforceIndependence = false (default here), the same runs are used at
% all time points, so each stage needs levels(k) runs per parameter set
% instead of levels(k)*numTimes.
%
% KDE bandwidths (bandwidthMode = 'local'): reference kernel widths are
% chosen once per stage, fitted time and species at the starting
% parameters (median of the full 'lcv' widths over bandwidthReps SSA
% repeats); every likelihood evaluation (and the calibration) then
% re-chooses the widths by cross-validation over 0.7-1.4 x the reference
% widths (5 candidates instead of 21). 'lcv' = full search every time;
% 'fixed' = reference widths without search (fastest but noisier).
%
% Outlier-robust likelihood (outlierFraction = eps > 0): each cell's
% SSA/KDE probability is replaced by (1-eps)*P_kde + eps*q, with q uniform
% over the observed range of nuclear and cytoplasmic counts. A few dozen
% cells far in the model's tails otherwise make up most of the noise of
% the SSA likelihood (see the per-cell diagnostic); the background bounds
% their influence. This assumes an outlier fraction eps in the data. The
% same eps is used in the transcription-site likelihood
% (computeTSlikelihood: P_obs = (1-eps)*P_model + eps*q, with q uniform over
% the observable TS counts 0 and TSthresh..50). A
% saved calibration or tuning state made with another eps is not reused.
%
% Calibration: the SSA likelihood noise sigma(N) is calibrated once, at
% the first starting point, and reused for all runs (optionally
% recalibrated every recalibrateEvery acceptances within a run).
%
% Parameter file (never Parameters_Final, which is only read):
%   parsFile ('savedParameters/Parameters_Final_DUSP1_w_GR_Binning') is
%   the working parameter file. If it does not exist, it is created with
%     * GRpars = GR parameters fitted to the GR intensities together with
%       the bin width (GRBinFits.<grBinModel>.pars(1:end-1) in grBinFile,
%       'savedParameters/Parameters_Final_GRBinning'); the fitted bin
%       width is stored in GRBinWidth;
%     * Pars_KON_KOFF_KR_zeroOFF, FinalModel_allButGR and
%       covPost_FinalModel_allButGR from sourceParsFile
%       ('savedParameters/Parameters_Final', read only), with the
%       GR-input parameters (24 kon_1, 25 mkoff, 26 kr_on_1) multiplied by
%       grInputScale (1.5) to account for the lower GR levels in units of
%       the fitted bin width (stored in GRInputScale).
%   The free parameters are read
%   from, and after each run written back to, the variables <parName> and
%   covPost_<parName>, where parName = 'FinalModel_allButGR' (joint data)
%   or 'FinalModel_allButGR_marginals' (marginalsOnly) -- the same naming
%   as GR_DUSP1_Fitting_Final, so the file can be used as its parameter
%   file. Writing to any file named Parameters_Final is an error.
%
% Output files (in outDir):
%   SeqABC_<tag>_calibration_<lossFun>.mat -- calibration
%   SeqABC_<tag>_tuning.mat   -- covPost, scaleMH, current parameters and
%                                tuning state (used to resume)
%   SeqABC_<tag>_chain.mat    -- progress of the current run (resumes an
%                                interrupted run; deleted when it ends)
%   SeqABC_<tag>_MHResults.mat -- final chain (MHResultsAll)
% Rerunning resumes from the saved tuning state and chain.
%
% New-best plots (plotNewBest, default true): whenever a run accepts a
% parameter set with a higher estimated log posterior than any earlier in
% that run, the data vs model marginals and means/stds vs time are
% plotted from the SSA runs already made for it and saved (overwritten)
% as <figDir>/SeqABC_<tag>_best_marginals.(fig/png) and
% <figDir>/SeqABC_<tag>_best_meanStd.(fig/png).
arguments
    opts.parsFile = 'savedParameters/Parameters_Final_DUSP1_w_GR_Binning'
    opts.sourceParsFile = 'savedParameters/Parameters_Final'
    opts.grBinFile = 'savedParameters/Parameters_Final_GRBinning'
    opts.grBinModel = 'combinedGRModel'
    opts.grInputScale (1,1) double {mustBePositive} = 1.5
    opts.modelLibrary = 'savedParameters/GRDusp1ModelLibrary'
    opts.dataFileDusp1 = 'Data/DUSP1_SSITcellresults_Final_Sep18.csv'
    opts.marginalsOnly (1,1) logical = false
    opts.kdeBiasCorrection {mustBeMember(opts.kdeBiasCorrection,{'none','jensen'})} = 'jensen'
    opts.outlierFraction (1,1) double = 1e-3
    opts.levels (1,:) double = [2500 10000 40000]
    opts.enforceIndependence (1,1) logical = false
    opts.bandwidthMode = 'local'
    opts.bandwidthReps (1,1) double = 3
    opts.z (1,1) double = 2.5
    opts.calibrationReps (1,1) double = 10
    opts.recalibrateEvery (1,1) double = Inf
    opts.includeTS (1,1) logical = true
    opts.usePrior (1,1) logical = true
    opts.tuneSamples (1,1) double = 200
    opts.tuneThin (1,1) double = 5
    opts.finalSamples (1,1) double = 10000
    opts.finalThin (1,1) double = 5
    opts.maxTuningRounds (1,1) double = 20
    opts.outDir = 'savedResults/SSALikelihood'
    opts.tag = ''
    opts.verbose (1,1) logical = true
    opts.plotNewBest (1,1) logical = true
    opts.figDir = 'savedFigures/SSALikelihood'
end
if isempty(opts.tag)
    opts.tag = 'TVandSat';
    if opts.marginalsOnly, opts.tag = 'TVandSat_marginals'; end
end
if ~isfolder(opts.outDir), mkdir(opts.outDir); end
lossFun = 'kde_loglik';
if opts.marginalsOnly, lossFun = 'kde_marginal_loglik'; end
calFile    = fullfile(opts.outDir, sprintf('SeqABC_%s_calibration_%s.mat', opts.tag, lossFun));
tuningFile = fullfile(opts.outDir, sprintf('SeqABC_%s_tuning.mat', opts.tag));
chainFile  = fullfile(opts.outDir, sprintf('SeqABC_%s_chain.mat', opts.tag));
resultFile = fullfile(opts.outDir, sprintf('SeqABC_%s_MHResults.mat', opts.tag));
figPrefix = '';
if opts.plotNewBest
    if ~isfolder(opts.figDir), mkdir(opts.figDir); end
    figPrefix = fullfile(opts.figDir, sprintf('SeqABC_%s_best', opts.tag));
end

% Working parameter file (created from sourceParsFile; never Parameters_Final).
parsFile = char(opts.parsFile);
if ~endsWith(parsFile, '.mat'), parsFile = [parsFile, '.mat']; end
assertNotParametersFinal(parsFile);
parName = 'FinalModel_allButGR';
if opts.marginalsOnly, parName = 'FinalModel_allButGR_marginals'; end
covName = ['covPost_', parName];
initWorkingParsFile(parsFile, opts.sourceParsFile, opts.grBinFile, opts.grBinModel, ...
    opts.grInputScale, parName, covName);

% Model, TS likelihood, and prior (as in GR_DUSP1_Fitting_Final, model 2).
[SSAModel, ModelTS, allData, indsFreePars, valsFreePars, log10PriorMean, log10PriorStd] = ...
    setupTVandSat(parsFile, opts.modelLibrary, opts.dataFileDusp1, parName);
TSLikelihood = @(pars) computeTSlikelihood(pars, ModelTS, allData, 50, true, opts.outlierFraction);
pm = log10PriorMean(indsFreePars); pm = pm(:);
ps = log10PriorStd(indsFreePars);  ps = ps(:);
logPriorDUSP1 = @(x) -sum((log10(x(:)) - pm).^2 ./ (2*ps.^2));
if opts.includeTS && opts.usePrior
    logExtra = @(x) safeEval(TSLikelihood, x) + logPriorDUSP1(x);
elseif opts.includeTS
    logExtra = @(x) safeEval(TSLikelihood, x);
elseif opts.usePrior
    logExtra = logPriorDUSP1;
else
    logExtra = @(x) 0;
end

% Run once to make sure files are up to date.
SSAModel=SSAModel.solve(solver='ssa');

% Tuning state: resume if saved, otherwise start as GR_DUSP1_Fitting_Final.
nPars = numel(valsFreePars);
if isfile(tuningFile)
    S = load(tuningFile);
    covPost = S.covPost; scaleMH = S.scaleMH; valsFreePars = S.valsFreePars;
    isTuning = S.isTuning; isFinal = S.isFinal; iRound = S.iRound;
    calibration = S.calibration;
    if ~isempty(calibration) && (~isequal(calibration.levels(:)', opts.levels(:)') || ...
            ~ssit.parest.calibrationOutlierMatches(calibration, opts.outlierFraction, []))
        fprintf('Saved calibration in %s does not match levels/outlierFraction; it will be redone.\n', tuningFile);
        calibration = [];
    end
    fprintf('Resuming sequential ABC tuning from %s (round %d, scaleMH = %g).\n', tuningFile, iRound, scaleMH);
else
    P = load(parsFile, covName);
    if isfield(P, covName)
        covPost = P.(covName);
    else
        covPost = 0.0001*eye(nPars);
    end
    scaleMH = 1; isTuning = true; isFinal = false; iRound = 0; calibration = [];
end

while isTuning || isFinal
    if ~isFinal
        nSamples = opts.tuneSamples; thin = opts.tuneThin;
    else
        nSamples = opts.finalSamples; thin = opts.finalThin;
    end
    covProp = scaleMH*covPost + 1e-6*eye(nPars);
    covProp = (covProp + covProp')/2;
    proposal = @(x) x + mvnrnd(0*x, covProp);
    if ~isfile(chainFile)
        % New run (otherwise an interrupted run is resumed from chainFile).
        iRound = iRound + 1;
        save(tuningFile, 'covPost', 'scaleMH', 'valsFreePars', 'isTuning', 'isFinal', ...
            'iRound', 'calibration', 'indsFreePars');
    end
    fprintf('Sequential ABC %s run %d: %d samples (thin %d), scaleMH = %g.\n', ...
        ternaryStr(isFinal,'final','tuning'), iRound, nSamples, thin, scaleMH);

    abcArgs = {'numberOfSamples', nSamples, 'thin', thin, 'burnIn', 0, ...
        'proposalDistribution', proposal, 'logPrior', logExtra, ...
        'marginalsOnly', opts.marginalsOnly, 'kdeBiasCorrection', opts.kdeBiasCorrection, ...
        'kdeOutlierFraction', opts.outlierFraction, ...
        'levels', opts.levels, 'enforceIndependence', opts.enforceIndependence, ...
        'bandwidthMode', opts.bandwidthMode, 'bandwidthReps', opts.bandwidthReps, ...
        'z', opts.z, 'calibrationReps', opts.calibrationReps, ...
        'recalibrateEvery', opts.recalibrateEvery, 'saveFile', chainFile, 'verbose', opts.verbose, ...
        'plotNewBest', opts.plotNewBest, 'newBestFigurePrefix', figPrefix};
    if isempty(calibration)
        abcArgs = [abcArgs, {'calibrationFile', calFile}];
    else
        abcArgs = [abcArgs, {'calibration', calibration}];
    end
    [parsBest, MHResultsAll] = SSAModel.runSequentialABC(valsFreePars, abcArgs{:});
    calibration = MHResultsAll.calibration;
    if isfile(chainFile), delete(chainFile); end

    % Update the proposal only if there are more accepted samples than
    % parameters (as in GR_DUSP1_Fitting_Final).
    acc = MHResultsAll.mhAcceptance;
    if acc*numel(MHResultsAll.mhValue) > size(MHResultsAll.mhSamples,2)
        covPost = cov(MHResultsAll.mhSamples);
        valsFreePars = parsBest;
        % Save to the working parameter file (as GR_DUSP1_Fitting_Final
        % saves covPost_<parName> and <parName> to its parameter file).
        toSave = struct(parName, valsFreePars, covName, covPost);
        assertNotParametersFinal(parsFile);
        save(parsFile, '-struct', 'toSave', '-append');
    end
    fprintf('Run %d acceptance rate: %.3f\n', iRound, acc);

    if isFinal
        save(resultFile, 'MHResultsAll', 'covPost', 'scaleMH', 'indsFreePars', 'opts');
        fprintf('Sequential ABC complete; results in %s.\n', resultFile);
        isTuning = false; isFinal = false;
    elseif acc < 0.05
        scaleMH = scaleMH*0.2;
    elseif acc < 0.2
        scaleMH = scaleMH*0.75;
    elseif acc > 0.6
        scaleMH = scaleMH*3;
    elseif acc > 0.45
        scaleMH = scaleMH*1.5;
    else
        isFinal = true;
    end
    if isTuning && ~isFinal && iRound >= opts.maxTuningRounds
        warning('Proposal tuning did not converge in %d runs; starting the final run.', iRound);
        isFinal = true;
    end
    save(tuningFile, 'covPost', 'scaleMH', 'valsFreePars', 'isTuning', 'isFinal', ...
        'iRound', 'calibration', 'indsFreePars');
end
out = struct('MHResultsAll', MHResultsAll, 'covPost', covPost, 'scaleMH', scaleMH, ...
    'pars', valsFreePars, 'indsFreePars', indsFreePars, 'calibration', calibration, ...
    'resultFile', resultFile, 'parsFile', parsFile);
end

function assertNotParametersFinal(file)
% Refuse to write to the original parameter file.
[~, name] = fileparts(char(file));
if strcmp(name, 'Parameters_Final')
    error('SSAbasedLogLikelihoodComparisons:protectedFile', ...
        'Writing to %s is not allowed; use a separate parameter file.', file);
end
end

function initWorkingParsFile(parsFile, sourceParsFile, grBinFile, grBinModel, grInputScale, parName, covName)
% Create the working parameter file: GR parameters from the GR binning fit
% (with the fitted bin width), DUSP1 parameters from the (read-only)
% source file with the GR-input parameters rescaled. Add <parName> /
% covPost_<parName> if they are missing. A file that is still a plain
% copy of the source (no GRBinWidth) is rebuilt; any other file without
% GRBinWidth is left alone and reported as an error.
assertNotParametersFinal(parsFile);
baseVars = {'GRpars','Pars_KON_KOFF_KR_zeroOFF','FinalModel_allButGR','covPost_FinalModel_allButGR'};
build = ~isfile(parsFile);
if ~build && ~ismember('GRBinWidth', who('-file', parsFile))
    W = load(parsFile);
    Src = load(sourceParsFile, baseVars{:});
    if isequal(sort(fieldnames(W)), sort(baseVars')) && ...
            all(cellfun(@(v) isequal(W.(v), Src.(v)), baseVars))
        build = true;
        fprintf('%s is a plain copy of %s; rebuilding it with the GR binning parameters.\n', ...
            parsFile, sourceParsFile);
    else
        error('SSAbasedLogLikelihoodComparisons:workingFile', ...
            ['%s exists but was not built from the GR binning fit (no GRBinWidth) ', ...
            'and is not a plain copy of %s. Move it or delete it to rebuild.'], parsFile, sourceParsFile);
    end
end
if build
    S = load(sourceParsFile, baseVars{:});
    B = load(grBinFile, 'GRBinFits');
    if ~isfield(B.GRBinFits, grBinModel)
        error('%s has no GRBinFits.%s; run GRBinning.m first.', grBinFile, grBinModel);
    end
    F = B.GRBinFits.(grBinModel);
    S.GRpars = F.pars(1:end-1);
    S.GRBinWidth = F.pars(end);
    S.GRBinModel = grBinModel;
    S.GRParsSource = grBinFile;
    % Rescale the GR-input parameters (24 kon_1, 25 mkoff, 26 kr_on_1).
    sets = {'Pars_KON_KOFF_KR_zeroOFF', [1,2,3,4,19,24,25,26]
            'FinalModel_allButGR',      [1,2,3,4,14,15,16,19,20,21,22,24,25,26]};
    for i = 1:size(sets,1)
        k = ismember(sets{i,2}, 24:26);
        S.(sets{i,1})(k) = grInputScale*S.(sets{i,1})(k);
    end
    S.GRInputScale = grInputScale;
    assertNotParametersFinal(parsFile);
    save(parsFile, '-struct', 'S');
    fprintf(['Created %s: GRpars from %s (GRBinFits.%s, bin width %.2f); DUSP1 parameters ', ...
        'from %s with GR-input parameters x %g.\n'], parsFile, grBinFile, grBinModel, ...
        S.GRBinWidth, sourceParsFile, grInputScale);
end
present = who('-file', parsFile);
missing = struct();
if ~ismember(parName, present)
    S = load(parsFile, 'FinalModel_allButGR');
    missing.(parName) = S.FinalModel_allButGR;
end
if ~ismember(covName, present) && ismember('covPost_FinalModel_allButGR', present)
    S = load(parsFile, 'covPost_FinalModel_allButGR');
    missing.(covName) = S.covPost_FinalModel_allButGR;
end
if ~isempty(fieldnames(missing))
    save(parsFile, '-struct', 'missing', '-append');
end
end

function [SSAModel, ModelTS, allData, indsFreePars, valsFreePars, log10PriorMean, log10PriorStd] = ...
    setupTVandSat(parsFile, modelLibrary, dataFileDusp1, parName)
% Model 2 (TVandSat) of GR_DUSP1_Fitting_Final: SSA model with nuclear and
% cytoplasmic DUSP1 data at 100 nM Dex, TS model and data, free parameters
% and prior, all set as in that function.
[SSAModel, log10PriorMean, log10PriorStd] = ...
    dusp1ModelLibrary_Final('fullSSAModel_100_2StateGen', false, modelLibrary);
DUSP1_Model_100nM = dusp1ModelLibrary_Final('ModelDUSP1_100nM', false, modelLibrary);
indsFreePars = [1,2,3,4,14,15,16,19,20,21,22,24,25,26];
zeroPars = 23;
P = load(parsFile, 'GRpars', 'Pars_KON_KOFF_KR_zeroOFF', parName);
valsFreePars = P.(parName)(:)';

SSAModel.parameters(5:12,2) = num2cell(P.GRpars);
SSAModel.parameters([1,2,3,4,19,24,25,26],2) = num2cell(P.Pars_KON_KOFF_KR_zeroOFF);
SSAModel.parameters(zeroPars,2) = num2cell(zeros(size(zeroPars)));
SSAModel.parameters(indsFreePars,2) = num2cell(valsFreePars);
SSAModel.fittingOptions.modelVarsToFit = indsFreePars;
SSAModel.solutionScheme = 'SSA';
if isunix && ~ismac
    % Same workaround as GR_DUSP1_Fitting_Final for nonlinear propensities.
    SSAModel.ssaOptions.useC = false;
end
SSAModel = SSAModel.solve(solver='ssa');   % compile the SSA once

ModelTS = DUSP1_Model_100nM;
ModelTS.dataSet = [];
ModelTS.parameters(5:12,2) = num2cell(P.GRpars);
ModelTS.parameters([1,2,3,4,19,24,25,26],2) = num2cell(P.Pars_KON_KOFF_KR_zeroOFF);
ModelTS.parameters(zeroPars,2) = num2cell(zeros(size(zeroPars)));
allData = readtable(dataFileDusp1);
allData = allData(allData.cyto_area_px>=12593 & allData.cyto_area_px<=17685 & isnan(allData.time_tpl),:);
ModelTS.tSpan = unique(allData.time)';
ModelTS.fittingOptions.modelVarsToFit = indsFreePars;
ModelTS.parameters(indsFreePars,2) = num2cell(valsFreePars);
end

function v = safeEval(fun, x)
% Evaluate a log-likelihood; failures (e.g., FSP errors) give -Inf (reject).
try
    v = fun(x);
    if ~isfinite(v), v = -Inf; end
catch err
    warning('SSAbasedLogLikelihoodComparisons:TSfailed', 'TS likelihood failed: %s', err.message);
    v = -Inf;
end
end

function s = ternaryStr(c, a, b)
if c, s = a; else, s = b; end
end
