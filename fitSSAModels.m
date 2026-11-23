function out = fitSSAModels(model, opts)
% fitSSAModels - Fit one of the full SSA DUSP1 models to the nuclear and
% cytoplasmic DUSP1 mRNA data (100 nM Dex) and the transcription-site
% data, using the sequential SSA-likelihood MH sampler
% (SSIT.runSequentialABC) with the MH settings and proposal tuning of
% GR_DUSP1_Fitting_Final. Adapted from seqABCTVandSat in
% SSAbasedLogLikelihoodComparisons.m.
%
%   out = fitSSAModels(model)
%   out = fitSSAModels(model, marginalsOnly=true, finalSamples=2000, ...)
%
% model -- number or name (as opts.FullModels in Main.m):
%   1 'TVandSat'    time-varying, saturating TTP degradation
%                   (GR_DUSP1_Fitting_Final model 2, FinalModel_allButGR)
%   2 'NoTV'        no time variation (model 4, FinalModelNoTimeVariation)
%   3 'NoSat'       no saturation (model 3, FinalModelNoSaturation)
%   4 'LinDecay'    first-order cytoplasmic decay (model 8,
%                   Pars_first_order_decay; no prior, as in that model)
%   5 'Mechanistic' extended mechanistic model with TTP (model 10,
%                   Pars_GR_DUSP1_TTP_Model_All)
% Free parameters, fixed-to-zero parameters, prior and starting values are
% those of the corresponding model in GR_DUSP1_Fitting_Final.
%
% Mechanistic model: the observed cytoplasmic count is rCyt + rCyt_complex
% (transformSSAsoln in GR_DUSP1_Fitting_Final). runSequentialABC compares
% model species directly, so a passive reporter species rCytTotal is
% added whose stoichiometry is the sum of the rows of rCyt and
% rCyt_complex (and whose initial value is their sum). It does not appear
% in any propensity, so the dynamics are unchanged and rCytTotal equals
% rCyt + rCyt_complex in every trajectory. The cytoplasmic data are linked
% to rCytTotal.
%
% Target (log posterior, natural log):
%   log L_SSA(nuclear/cytoplasmic data)  -- SSA/KDE estimate, sequential test
% + log L_TS(transcription-site data)    -- FSP (computeTSlikelihood)
% + log prior (if used)                  -- as in computeTotalObjFun:
%       -sum((log10(x)-log10PriorMean).^2 ./ (2*log10PriorStd.^2))
% The TS likelihood and prior are deterministic, so they enter the MH
% threshold (runSequentialABC's logPrior input) and only the SSA term is
% estimated sequentially. Unlike GR_DUSP1_Fitting_Final, the terms are not
% weighted (weights_ABC there); all are log-likelihoods/log-priors.
%
% Outlier-robust likelihood (outlierFraction = eps > 0): each cell's
% probability is (1-eps)*P_model + eps*q, with q uniform over the observed
% range, in both the SSA/KDE likelihood (nuclear and cytoplasmic counts)
% and the TS likelihood (observable TS counts 0 and 4..50).
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
% SSA settings: stages of levels runs per parameter set; with
% enforceIndependence = false the same runs are used at all time points.
% KDE bandwidths (bandwidthMode = 'local') are chosen near reference
% widths found once at calibration. The SSA likelihood noise sigma(N) is
% calibrated once per model at the first starting point.
%
% Parameter file (never Parameters_Final, which is only read):
%   parsFile ('savedParameters/Parameters_Final_DUSP1_w_GR_Binning') is
%   shared by all models. If it does not exist, it is created with GRpars
%   from the GR binning fit (grBinFile, with the fitted bin width in
%   GRBinWidth) and the DUSP1 parameters from sourceParsFile, with the
%   GR-input parameters (24 kon_1, 25 mkoff, 26 kr_on_1) multiplied by
%   grInputScale. The free parameters of each model are read from, and
%   after each run written back to, <parName> and covPost_<parName>
%   (parName as listed above, with '_marginals' appended if
%   marginalsOnly). If <parName> is missing, it is copied from
%   sourceParsFile with the same GR-input rescaling.
%
% Output files (in outDir, tag = model name [+ '_marginals']):
%   SeqABC_<tag>_calibration_<lossFun>.mat, SeqABC_<tag>_tuning.mat,
%   SeqABC_<tag>_chain.mat (resume file), SeqABC_<tag>_MHResults.mat.
% New-best plots: <figDir>/SeqABC_<tag>_best_marginals / _meanStd.
arguments
    model
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
    opts.usePrior = []          % [] = as in GR_DUSP1_Fitting_Final for this model
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
cfg = modelConfig(model);
if isempty(opts.usePrior), opts.usePrior = cfg.usePrior; end
if isempty(opts.tag)
    opts.tag = cfg.name;
    if opts.marginalsOnly, opts.tag = [cfg.name, '_marginals']; end
end
fprintf('Fitting SSA model %d (%s): %d free parameters, prior %s.\n', ...
    cfg.number, cfg.name, numel(cfg.indsFreePars), ternaryStr(opts.usePrior,'on','off'));

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
parName = cfg.parName;
if opts.marginalsOnly, parName = [cfg.parName, '_marginals']; end
covName = ['covPost_', parName];
initWorkingParsFile(parsFile, opts.sourceParsFile, opts.grBinFile, opts.grBinModel, ...
    opts.grInputScale, cfg, parName, covName);

% Model, TS likelihood, and prior (as in GR_DUSP1_Fitting_Final).
[SSAModel, ModelTS, allData, valsFreePars, log10PriorMean, log10PriorStd] = ...
    setupModel(cfg, parsFile, opts.modelLibrary, opts.dataFileDusp1, parName);
indsFreePars = cfg.indsFreePars;
TSLikelihood = @(pars) computeTSlikelihood(pars, ModelTS, allData, 50, true, opts.outlierFraction);
if opts.usePrior
    if isempty(log10PriorMean)
        error('fitSSAModels:noPrior', 'The model library returned no prior for %s.', cfg.libraryName);
    end
    pm = log10PriorMean(indsFreePars); pm = pm(:);
    ps = log10PriorStd(indsFreePars);  ps = ps(:);
    logPriorDUSP1 = @(x) -sum((log10(x(:)) - pm).^2 ./ (2*ps.^2));
end
if opts.includeTS && opts.usePrior
    logExtra = @(x) safeEval(TSLikelihood, x) + logPriorDUSP1(x);
elseif opts.includeTS
    logExtra = @(x) safeEval(TSLikelihood, x);
elseif opts.usePrior
    logExtra = logPriorDUSP1;
else
    logExtra = @(x) 0;
end

% Tuning state: resume if saved, otherwise start as GR_DUSP1_Fitting_Final.
nPars = numel(valsFreePars);
if isfile(tuningFile)
    S = load(tuningFile);
    if ~isequal(S.indsFreePars(:)', indsFreePars(:)')
        error('fitSSAModels:tuningFile', ...
            '%s was saved for other free parameters; move or delete it.', tuningFile);
    end
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
    fprintf('Sequential ABC (%s) %s run %d: %d samples (thin %d), scaleMH = %g.\n', ...
        cfg.name, ternaryStr(isFinal,'final','tuning'), iRound, nSamples, thin, scaleMH);

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
        abcArgs = [abcArgs, {'calibrationFile', calFile}]; %#ok<AGROW>
    else
        abcArgs = [abcArgs, {'calibration', calibration}]; %#ok<AGROW>
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
        toSave = struct(parName, valsFreePars, covName, covPost);
        assertNotParametersFinal(parsFile);
        save(parsFile, '-struct', 'toSave', '-append');
    end
    fprintf('Run %d acceptance rate: %.3f\n', iRound, acc);

    if isFinal
        save(resultFile, 'MHResultsAll', 'covPost', 'scaleMH', 'indsFreePars', 'opts', 'cfg');
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
out = struct('model', cfg.name, 'MHResultsAll', MHResultsAll, 'covPost', covPost, ...
    'scaleMH', scaleMH, 'pars', valsFreePars, 'indsFreePars', indsFreePars, ...
    'parName', parName, 'calibration', calibration, 'resultFile', resultFile, ...
    'parsFile', parsFile);
end

% =========================================================================
function cfg = modelConfig(model)
% Settings of each SSA model, as in GR_DUSP1_Fitting_Final.
names = {'TVandSat','NoTV','NoSat','LinDecay','Mechanistic'};
if ischar(model) || isstring(model)
    k = find(strcmpi(names, char(model)));
    if isempty(k)
        error('fitSSAModels:model', 'Unknown model ''%s''. Use one of: %s.', ...
            char(model), strjoin(names, ', '));
    end
else
    k = model;
    if ~isscalar(k) || ~ismember(k, 1:numel(names))
        error('fitSSAModels:model', 'model must be 1-%d or one of: %s.', ...
            numel(names), strjoin(names, ', '));
    end
end
cfg = struct('number', k, 'name', names{k}, 'libraryName', '', ...
    'indsFreePars', [], 'zeroPars', [], 'parName', '', 'usePrior', true, ...
    'linearDecay', false, 'addCytReporter', false, 'tSpanFromMinus500', false, ...
    'fittingFinalModel', NaN);
switch k
    case 1 % TVandSat (GR_DUSP1_Fitting_Final model 2)
        cfg.libraryName = 'fullSSAModel_100_2StateGen';
        cfg.indsFreePars = [1,2,3,4,14,15,16,19,20,21,22,24,25,26];
        cfg.zeroPars = 23;
        cfg.parName = 'FinalModel_allButGR';
        cfg.fittingFinalModel = 2;
    case 2 % NoTV (model 4)
        cfg.libraryName = 'altSSAModel_NoTimeVariation_100';
        cfg.indsFreePars = [1,2,3,4,14,15,16,19,20,21,24,25,26];
        cfg.zeroPars = 23;
        cfg.parName = 'FinalModelNoTimeVariation';
        cfg.fittingFinalModel = 4;
    case 3 % NoSat (model 3)
        cfg.libraryName = 'altSSAModel_NoSaturation_100';
        cfg.indsFreePars = [1,2,3,4,14,15,19,20,21,22,24,25,26];
        cfg.zeroPars = 23;
        cfg.parName = 'FinalModelNoSaturation';
        cfg.fittingFinalModel = 3;
    case 4 % LinDecay (model 8): first-order decay, no prior
        cfg.libraryName = 'fullSSAModel_100_2StateGen';
        cfg.indsFreePars = [1,2,3,4,14,19,24,25,26];
        cfg.zeroPars = [15,16,21,22,23];
        cfg.parName = 'Pars_first_order_decay';
        cfg.usePrior = false;
        cfg.linearDecay = true;
        cfg.fittingFinalModel = 8;
    case 5 % Mechanistic (model 10): TTP model with prior
        cfg.libraryName = 'GR_DUSP1_TTP_Model_2';
        cfg.indsFreePars = [1,2,3,4,14,15,19,24,25,26,27:36];
        cfg.zeroPars = [20,21,22,23];
        cfg.parName = 'Pars_GR_DUSP1_TTP_Model_All';
        cfg.addCytReporter = true;
        cfg.tSpanFromMinus500 = true;
        cfg.fittingFinalModel = 10;
end
end

function assertNotParametersFinal(file)
% Refuse to write to the original parameter file.
[~, name] = fileparts(char(file));
if strcmp(name, 'Parameters_Final')
    error('fitSSAModels:protectedFile', ...
        'Writing to %s is not allowed; use a separate parameter file.', file);
end
end

function initWorkingParsFile(parsFile, sourceParsFile, grBinFile, grBinModel, grInputScale, cfg, parName, covName)
% Create the working parameter file (as in SSAbasedLogLikelihoodComparisons):
% GR parameters from the GR binning fit, DUSP1 parameters from the
% (read-only) source file with the GR-input parameters rescaled. Then add
% this model's <parName> / covPost_<parName> if they are missing.
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
        error('fitSSAModels:workingFile', ...
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

% This model's free parameters (and proposal covariance) in the working
% file: copied from the source file (GR-input parameters rescaled), or,
% for the '_marginals' variant, from the joint-data values if present.
present = who('-file', parsFile);
missing = struct();
if ~ismember(parName, present)
    if ~strcmp(parName, cfg.parName) && ismember(cfg.parName, present)
        S = load(parsFile, cfg.parName);
        missing.(parName) = S.(cfg.parName);
    else
        S = load(sourceParsFile, cfg.parName);
        v = S.(cfg.parName);
        if numel(v) ~= numel(cfg.indsFreePars)
            error('fitSSAModels:parSize', '%s in %s has %d entries, but %s has %d free parameters.', ...
                cfg.parName, sourceParsFile, numel(v), cfg.name, numel(cfg.indsFreePars));
        end
        k = ismember(cfg.indsFreePars, 24:26);
        v(k) = grInputScale*v(k);
        missing.(parName) = v;
        fprintf('Added %s to %s (from %s, GR-input parameters x %g).\n', ...
            parName, parsFile, sourceParsFile, grInputScale);
    end
end
if ~ismember(covName, present)
    srcCov = ['covPost_', cfg.parName];
    if ismember(srcCov, present)
        S = load(parsFile, srcCov);
        missing.(covName) = S.(srcCov);
    elseif ismember(srcCov, who('-file', sourceParsFile))
        S = load(sourceParsFile, srcCov);
        missing.(covName) = S.(srcCov);
    end
end
if ~isempty(fieldnames(missing))
    assertNotParametersFinal(parsFile);
    save(parsFile, '-struct', 'missing', '-append');
end
end

function [SSAModel, ModelTS, allData, valsFreePars, log10PriorMean, log10PriorStd] = ...
    setupModel(cfg, parsFile, modelLibrary, dataFileDusp1, parName)
% SSA model with nuclear and cytoplasmic DUSP1 data at 100 nM Dex, TS
% model and data, free parameters and prior, set as in
% GR_DUSP1_Fitting_Final for model cfg.fittingFinalModel.
[SSAModel, log10PriorMean, log10PriorStd] = ...
    dusp1ModelLibrary_Final(cfg.libraryName, false, modelLibrary);
DUSP1_Model_100nM = dusp1ModelLibrary_Final('ModelDUSP1_100nM', false, modelLibrary);
indsFreePars = cfg.indsFreePars;
zeroPars = cfg.zeroPars;
P = load(parsFile, 'GRpars', 'Pars_KON_KOFF_KR_zeroOFF', parName);
valsFreePars = P.(parName)(:)';

if cfg.linearDecay
    % First-order cytoplasmic decay (GR_DUSP1_Fitting_Final model 8).
    SSAModel.propensityFunctions(10) = {'degCyt0*rCyt'};
end
if cfg.addCytReporter
    % Observed cytoplasmic count = rCyt + rCyt_complex: add a passive
    % reporter species rCytTotal (not used in any propensity) that tracks
    % this sum exactly, and link the cytoplasmic data to it.
    iCyt = find(strcmp(SSAModel.species, 'rCyt'));
    iCmplx = find(strcmp(SSAModel.species, 'rCyt_complex'));
    if isempty(iCyt) || isempty(iCmplx)
        error('fitSSAModels:species', '%s has no rCyt or rCyt_complex species.', cfg.libraryName);
    end
    SSAModel.species{end+1} = 'rCytTotal';
    SSAModel.initialCondition(end+1) = SSAModel.initialCondition(iCyt) + SSAModel.initialCondition(iCmplx);
    SSAModel.stoichiometry(end+1,:) = SSAModel.stoichiometry(iCyt,:) + SSAModel.stoichiometry(iCmplx,:);
    SSAModel = SSAModel.loadData(dataFileDusp1, ...
        {'rna','num_nuc_noTS'; 'rCytTotal','num_cyto_spots'}, ...
        {[],[],['(TAB.dex_conc==100|(TAB.dex_conc==0&TAB.time==0))', ...
        '&TAB.cyto_area_px>=12593&TAB.cyto_area_px<=17685&isnan(TAB.time_tpl)']}, {'replica'});
end
if cfg.tSpanFromMinus500
    SSAModel.tSpan = [-500, SSAModel.dataSet.times];
end

SSAModel.parameters(5:12,2) = num2cell(P.GRpars);
SSAModel.parameters([1,2,3,4,19,24,25,26],2) = num2cell(P.Pars_KON_KOFF_KR_zeroOFF);
SSAModel.parameters(zeroPars,2) = num2cell(zeros(size(zeroPars)));
SSAModel.parameters(indsFreePars,2) = num2cell(valsFreePars);
SSAModel.fittingOptions.modelVarsToFit = indsFreePars;
SSAModel.solutionScheme = 'SSA';
if isunix && ~ismac && cfg.fittingFinalModel < 10
    % Same workaround as GR_DUSP1_Fitting_Final for nonlinear propensities.
    SSAModel.ssaOptions.useC = false;
end
% Write and compile the SSA code for this (possibly modified) model once.
SSAModel.ssaOptions.computeFile = '';
SSAModel = SSAModel.solve(solver='ssa');

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
    warning('fitSSAModels:TSfailed', 'TS likelihood failed: %s', err.message);
    v = -Inf;
end
end

function s = ternaryStr(c, a, b)
if c, s = a; else, s = b; end
end
