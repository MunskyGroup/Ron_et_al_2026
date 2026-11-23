%% GRBinning.m
% Version of the GR model that predicts the measured (background-corrected)
% GR fluorescence intensities directly, using a binning PDO that maps the
% model's discrete GR "count" (bin index) to a uniform distribution of
% intensities within that bin.
%
% PDO definition (for the cytoplasm; the nucleus uses the same rules with
% the threshold and width divided by n2cRatio):
%   true count x = 0   ->  intensity uniform on [0, T)
%   true count x >= 1  ->  intensity uniform on [T+(x-1)*w, T+x*w)
% with T = 1087 and w = binWidth (a fitting variable).
%
% The current (discretized) approach corresponds to binning the data with
% the same rule (normGRcyt/normGRnuc columns). So the likelihood of the
% intensity data under this PDO equals the likelihood of the binned data
% minus log(bin width) for each cell, up to the intensity quantization below.
%
% NOTE on intensity quantization: SSIT stores data in a sparse tensor indexed
% by integer values, and the PDO-distorted joint (cyt x nuc) distribution
% has one entry per pair of observable intensity values. Intensities are
% therefore rounded to units of 'intensityUnit' (a.u.). The default of 10 a.u.
% is much finer than the bin widths (~82 and ~149 a.u.) and keeps the joint
% observed distribution manageable (~600 x 1000 per time point instead of
% ~6000 x 10000 at 1 a.u.).

clear
addpath(genpath('src'));

modelLibrary  = 'savedParameters/GRDusp1ModelLibrary';
savedParsFile = 'savedParameters/Parameters_Final';
dataFileGR    = 'Data/GR_SSITcellresults_Final_Sep18.csv';

%% Binning / PDO settings
binOpts.cytThreshold  = 1087;    % Cyt intensities below this map to bin 0.
binOpts.n2cRatio      = 0.5514;  % Nucleus-to-cytoplasm area ratio.
binOpts.intensityUnit = 10;      % Intensity quantization (a.u.) for data and PDO.
binOpts.dataColumns   = {'cytGR','CytGRCorrected';'nucGR','NucGRCorrected'};

%% GR model versions
% Library models, the labels used for them in Main.m, and their parameters
% fitted to the discretized data in savedParsFile.
GRModelNames  = {'combinedGRModel','combinedGRModel_CytDegOnly', ...
    'combinedGRModel_NucDegOnly','combinedGRModel_DexExport'};
% (combinedGRModel_CytDegOnly has no nuclear degradation, gg2 = 0;
% combinedGRModel_NucDegOnly has no cytoplasmic degradation, gg1 = 0.)
GRModelLabels = {'DexControlsGRImport_NucCytDeg','DexControlsGRImport_CytDeg', ...
    'DexControlsGRImport_NucDeg','DexControlsGRExport_NucCytDeg'};
savedGRParSet = {'GRpars','Pars_CytDegOnly','Pars_NucDegOnly','Pars_DexExport'};

%% Update the 100, 10, and 1 nM Dex models to include a PDO where:
% Cytoplasmic intensities below 1087 map to bin 0.
% The nucleus to cytoplasmic area ratio is n2cRatio=0.5514
% Nuclear intensities below 1087/n2cRatio map to bin 0.
% Cytoplasmic bins have a width 'binWidth' which is a fitting variable.
% Nuclear bins have a width 'binWidth/n2cRatio'.
% PDO maps within each bin are uniform over the bin widths.
% Data loaded into each model is from the columns "CytGRCorrected" and
% "NucGRCorrected".
% This uses the SSIT PDO type 'Intensity Binning - Parametrized', with
% binWidth as its (only) PDO parameter. See buildGRBinModel below.
binWidth = 82;   % Starting bin width.

%% Make an example plot of the results
% Start with a binWidth = 82.
% Load the best fit model from previous analysis that fit GR to the
% discretized data.
% Make a plot of the distorted data and model prediction using the PDO.
[combinedGRBinModel, combinedGRModel] = buildGRBinModel('combinedGRModel', ...
    binWidth, binOpts, modelLibrary, dataFileGR);
load(savedParsFile,'GRpars');
combinedGRModel    = combinedGRModel.updateModels(GRpars, false);
combinedGRBinModel = combinedGRBinModel.updateModels([GRpars, binWidth], false);

logLBinned    = combinedGRModel.computeTotalLogLikelihood(GRpars);
logLIntensity = combinedGRBinModel.computeTotalLogLikelihood([GRpars, binWidth]);
fprintf('Log-posterior, discretized GR data (normGR columns):        %g\n', logLBinned);
fprintf('Log-posterior, intensity data with binning PDO (w = %g):  %g\n', binWidth, logLIntensity);

figHandles = plotGRIntensityFits(combinedGRBinModel, binOpts);

%% Fit the GR Model and the Bidwidth Simultaneously.
% Fits each GR model version (GR parameters and binWidth together) with
% fminsearch in log space, in nFitRounds rounds. In each round, every model
% in modelsToFit gets maxIterPerRound iterations, restarting fminsearch
% from its best parameters so far. After each model's round, results are
% saved to savedParsFileBin (a copy of savedParsFile) in the struct
% GRBinFits, with one field per library model name:
%   label    - model label used in Main.m
%   parNames - names of the free parameters (GR parameters, then binWidth)
%   pars     - best parameters
%   logPost  - log posterior at pars
%   fitTrace - every evaluation: [round, eval #, seconds, -log posterior]
% and that model's figures are updated and saved to figDirBin. Models with
% results in GRBinFits continue from them; others start from their
% discretized-data fit in savedParsFile with binWidth = 82.
runFits          = true;
modelsToFit      = 1:4;     % indices into GRModelNames
nFitRounds       = 5;
maxIterPerRound  = 100;
savedParsFileBin = 'savedParameters/Parameters_Final_GRBinning';
figDirBin        = 'savedFigures/GRBinning';

if runFits
    if ~exist([savedParsFileBin,'.mat'],'file')
        copyfile([savedParsFile,'.mat'],[savedParsFileBin,'.mat']);
    end
    if ~exist(figDirBin,'dir')
        mkdir(figDirBin);
    end
    S = load(savedParsFileBin);
    if isfield(S,'GRBinFits')
        GRBinFits = S.GRBinFits;
    else
        GRBinFits = struct();
    end
    % Move results saved by an earlier version of this script (base model
    % only, as separate GRBin* variables) into GRBinFits.
    legacyVars = {'GRBinPars','GRBinParNames','GRBinLogPost','GRBinFitTrace'};
    if isfield(S,'GRBinPars')
        if ~isfield(GRBinFits,'combinedGRModel')
            GRBinFits.combinedGRModel = struct('label',GRModelLabels{1}, ...
                'parNames',{S.GRBinParNames}, 'pars',S.GRBinPars, ...
                'logPost',S.GRBinLogPost, 'fitTrace',S.GRBinFitTrace);
        end
        S = rmfield(S, legacyVars(isfield(S,legacyVars)));
        S.GRBinFits = GRBinFits;
        save(savedParsFileBin, '-struct', 'S');
    end
    clear S

    % Build binned versions of all models to fit and set starting values.
    binModels = cell(1,length(GRModelNames));
    for iM = modelsToFit
        name = GRModelNames{iM};
        binModels{iM} = buildGRBinModel(name, binWidth, binOpts, modelLibrary, dataFileGR);
        if ~isfield(GRBinFits,name)
            P = load(savedParsFile, savedGRParSet{iM});
            M1 = binModels{iM}.SSITModels{1};
            F.label    = GRModelLabels{iM};
            F.parNames = [M1.parameters(M1.fittingOptions.modelVarsToFit,1); {'binWidth'}]';
            F.pars     = [P.(savedGRParSet{iM}), binWidth];
            F.logPost  = NaN;
            F.fitTrace = zeros(0,4);
            GRBinFits.(name) = F;
        end
    end

    figHandlesFit = cell(1,length(GRModelNames));
    figConvergence = [];
    for iRound = 1:nFitRounds
        for iM = modelsToFit
            name = GRModelNames{iM};
            F = GRBinFits.(name);
            F.label = GRModelLabels{iM};
            if isempty(F.fitTrace)
                modelRound = 1; evalOffset = 0; timeOffset = 0;
            else
                modelRound = F.fitTrace(end,1) + 1;
                evalOffset = F.fitTrace(end,2); timeOffset = F.fitTrace(end,3);
            end
            fprintf('\n===== %s (%s), round %d: fminsearch, %d iterations =====\n', ...
                name, F.label, modelRound, maxIterPerRound);
            parsStart = F.pars;
            [F.pars, trace] = fitGRBinModel(binModels{iM}, log(F.pars), maxIterPerRound);
            F.logPost = -min(trace(:,3));
            F.fitTrace = [F.fitTrace; modelRound*ones(size(trace,1),1), ...
                trace(:,1)+evalOffset, trace(:,2)+timeOffset, trace(:,3)];
            GRBinFits.(name) = F;

            % Save parameters.
            save(savedParsFileBin, 'GRBinFits', '-append');
            fprintf('%s round %d: log posterior = %.1f, binWidth = %.2f (saved to %s)\n', ...
                name, modelRound, F.logPost, F.pars(end), savedParsFileBin);
            disp(array2table([parsStart(:), F.pars(:)], 'VariableNames', {'start','fit'}, ...
                'RowNames', F.parNames))

            % Update and save figures.
            binModels{iM} = binModels{iM}.updateModels(F.pars, false);
            figHandlesFit{iM} = plotGRIntensityFits(binModels{iM}, binOpts, figHandlesFit{iM}, name);
            for iFig = 1:length(figHandlesFit{iM})
                dex = binModels{iM}.SSITModels{iFig}.parameters{13,2};
                saveFigure(figHandlesFit{iM}{iFig}, ...
                    fullfile(figDirBin, sprintf('%s_GRIntensityFits_Dex%g', name, dex)));
            end
            figConvergence = plotFitConvergence(GRBinFits, figConvergence);
            saveFigure(figConvergence, fullfile(figDirBin, 'GRBinFitConvergence'));
            drawnow
        end
    end

    % Summary of all fitted models.
    names = fieldnames(GRBinFits);
    summary = table(cellfun(@(n)GRBinFits.(n).label, names, 'UniformOutput', false), ...
        cellfun(@(n)numel(GRBinFits.(n).pars), names), ...
        cellfun(@(n)GRBinFits.(n).logPost, names), ...
        cellfun(@(n)GRBinFits.(n).pars(end), names), ...
        cellfun(@(n)size(GRBinFits.(n).fitTrace,1), names), ...
        'VariableNames', {'label','nPars','logPost','binWidth','nEvals'}, 'RowNames', names);
    disp(summary)
end

%% ************************************************************************
%  Local functions
%  ************************************************************************
function [binSMM, SMM] = buildGRBinModel(modelName, binWidth, binOpts, modelLibrary, dataFileGR)
% Load GR model version 'modelName' from the library and return a version
% (binSMM) that fits the measured intensities through the binning PDO. SMM
% is the original model fit to discretized data. Each of the three Dex
% models shares the GR parameters and the bin width, which is the last
% free parameter (each model's free parameters are its model parameters
% followed by its PDO parameters). The prior applies to the GR parameters
% only.
SMM = dusp1ModelLibrary_Final(modelName,false,modelLibrary);
nMods = length(SMM.SSITModels);
binModels = cell(1,nMods);
for iMod = 1:nMods
    Model = SMM.SSITModels{iMod};
    % Reload the intensity data, using the same filtering conditions (Dex
    % concentration, excluded time points, etc.) as the discretized model.
    linkedSpecies = binOpts.dataColumns;
    for iS = 1:size(linkedSpecies,1)
        linkedSpecies(iS,2:3) = {'', ...
            sprintf('round(TAB.%s/%g)',linkedSpecies{iS,2},binOpts.intensityUnit)};
    end
    Model = Model.loadData(dataFileGR, linkedSpecies, Model.dataSet.conditions);
    Model.tSpan = Model.dataSet.times;
    binModels{iMod} = addGRIntensityPDO(Model, binWidth, binOpts);
end
nGRPars = length(SMM.parameters);
priorGR = SMM.parameterConstraints;
binSMM = SSITMultiModel(binModels, ...
    cellfun(@(p)[p,nGRPars+1], SMM.parameterIndices, 'UniformOutput', false), ...
    @(x)priorGR(x(1:nGRPars)), SMM.fspStateSpaces);
end

function Model = addGRIntensityPDO(Model, binWidth, binOpts)
% Attach the SSIT 'Intensity Binning - Parametrized' PDO to a GR model, with
% the cytoplasmic bin width as the single fittable PDO parameter C(1).
% Thresholds and widths are in units of binOpts.intensityUnit.
u = binOpts.intensityUnit;
r = binOpts.n2cRatio;
Model.pdoOptions.type = 'Intensity Binning - Parametrized';
Model.pdoOptions.props = struct('ParameterGuess', binWidth);
for iS = 1:length(Model.species)
    sp = num2str(iS);
    switch Model.species{iS}
        case 'cytGR'
            Model.pdoOptions.props.(['BinThresholdS',sp]) = binOpts.cytThreshold/u;
            Model.pdoOptions.props.(['BinWidthS',sp]) = @(C)C(1)/u;
        case 'nucGR'
            Model.pdoOptions.props.(['BinThresholdS',sp]) = binOpts.cytThreshold/r/u;
            Model.pdoOptions.props.(['BinWidthS',sp]) = @(C)C(1)/r/u;
        otherwise
            error('No intensity binning defined for species %s.', Model.species{iS})
    end
end
Model.fittingOptions.pdoVarsToFit = 1;
end

function [pars, trace] = fitGRBinModel(SMM, x0, maxIter)
% Fit all free parameters of the multimodel (GR parameters and binWidth)
% with fminsearch in log space for maxIter iterations. Returns the best
% parameters found and a trace of every evaluation:
% [evaluation #, seconds since start, -log posterior].
recordEval('reset');
objFun = @(x)recordEval(-SMM.computeTotalLogLikelihood(exp(x(:)')), x);
opts = optimset('Display','iter','MaxIter',maxIter,'MaxFunEvals',inf);
fminsearch(objFun, x0(:), opts);
trace = recordEval('get');
[~,iBest] = min(trace(:,3));
pars = exp(trace(iBest,4:end));
trace = trace(:,1:3);
end

function out = recordEval(f, x)
% Records every objective evaluation: [eval #, seconds, f, x'].
% recordEval('reset') clears, recordEval('get') returns the record, and
% recordEval(f,x) appends an evaluation and returns f.
persistent T tStart
if ischar(f)
    switch f
        case 'reset'
            T = []; tStart = tic;
        case 'get'
            out = T;
    end
    return
end
T(end+1,:) = [size(T,1)+1, toc(tStart), f, x(:)'];
out = f;
end

function fig = plotFitConvergence(GRBinFits, fig)
% Best -log posterior so far versus likelihood evaluations for each fitted
% model. Pass an existing figure handle to redraw it in place.
if isempty(fig) || ~isgraphics(fig)
    fig = figure('Name','GR + binWidth fit convergence','Position',[100 100 700 450]);
else
    clf(fig);
end
ax = axes(fig); hold(ax,'on');
names = fieldnames(GRBinFits);
leg = {};
for i = 1:length(names)
    T = GRBinFits.(names{i}).fitTrace;
    if ~isempty(T)
        stairs(ax, T(:,2), cummin(T(:,4)), 'LineWidth', 2);
        leg{end+1} = sprintf('%s (%.1f)', GRBinFits.(names{i}).label, min(T(:,4))); %#ok<AGROW>
    end
end
hold(ax,'off'); grid(ax,'on');
xlabel(ax,'likelihood evaluations');
ylabel(ax,'best -log posterior');
title(ax,'Joint fits of GR parameters and binWidth');
legend(ax, leg, 'Location','northeast', 'Interpreter','none');
end

function saveFigure(fig, baseName)
% Save a figure as .fig and .png.
savefig(fig, [baseName,'.fig']);
exportgraphics(fig, [baseName,'.png'], 'Resolution', 150);
end

function [t,v] = speciesBinGeometry(speciesName, binWidth, binOpts)
% Threshold t and bin width v (in units of binOpts.intensityUnit) for a
% species, matching addGRIntensityPDO (used to draw histogram bin edges).
switch speciesName
    case 'cytGR'
        scale = 1;
    case 'nucGR'
        scale = 1/binOpts.n2cRatio;
    otherwise
        error('No intensity binning defined for species %s.', speciesName)
end
t = scale*binOpts.cytThreshold/binOpts.intensityUnit;
v = scale*binWidth/binOpts.intensityUnit;
end

function figHandles = plotGRIntensityFits(combinedModel, binOpts, figHandles, modelName)
% For each Dex concentration, plot the measured cytoplasmic and nuclear GR
% intensity distributions (histograms with the PDO bin edges) against the
% PDO-distorted model predictions. Pass the figHandles from a previous call
% to redraw those figures in place.
arguments
    combinedModel
    binOpts
    figHandles = []
    modelName = ''
end
u = binOpts.intensityUnit;
nMods = length(combinedModel.SSITModels);
if isempty(figHandles)
    figHandles = cell(1,nMods);
end

% Common x-limits for all panels, spanning all intensity data (all Dex
% concentrations, times, and compartments), rounded out to 500 a.u.
allData = cellfun(@(m) reshape(cell2mat(m.dataSet.DATA(:,2:end)),[],1), ...
    combinedModel.SSITModels, 'UniformOutput', false);
allData = vertcat(allData{:})*u;
xLim = [500*floor(min(allData)/500), 500*ceil(max(allData)/500)];

for iMod = 1:nMods
    Model = combinedModel.SSITModels{iMod};
    Model = Model.solve(solver='fsp');
    [PDO,Model] = Model.generatePDO([],[],Model.Solutions.fsp);
    binWidth = Model.pdoOptions.props.ParameterGuess(1);
    dex = Model.parameters{strcmp(Model.parameters(:,1),'Dex0'),2};

    DATA = cell2mat(Model.dataSet.DATA);   % columns: time, species (model order)
    times = Model.dataSet.times;
    nT = length(times);
    nS = length(Model.species);

    if isempty(figHandles{iMod}) || ~isgraphics(figHandles{iMod})
        figHandles{iMod} = figure('Name',sprintf('%s GR intensity fits, %g nM Dex',modelName,dex), ...
            'Position',[50 50 230*nT 420]);
    else
        clf(figHandles{iMod});
    end
    tl = tiledlayout(figHandles{iMod}, nS, nT, 'TileSpacing','compact', 'Padding','compact');
    title(tl, sprintf('%s GR intensities, %g nM Dex (cyt bin width = %.1f)', modelName, dex, binWidth), ...
        'Interpreter','none');

    for iT = 1:nT
        [~,jT] = min(abs(Model.tSpan - times(iT)));
        P = double(Model.Solutions.fsp{jT}.p.data);
        for iS = 1:nS
            % Marginal of true counts and its PDO distortion. The PDO is a
            % tensor product, so marginals transform independently.
            pTrue = sum(P, setdiff(1:ndims(P), iS));
            pTrue = pTrue(:);
            C = PDO.conditionalPmfs{iS};
            nX = min(length(pTrue), size(C,2));
            pObs = C(:,1:nX)*pTrue(1:nX);

            [t,v] = speciesBinGeometry(Model.species{iS}, binWidth, binOpts);
            xData = DATA(DATA(:,1)==times(iT), iS+1)*u;
            edges = u*[0, t + v*(0:ceil((xLim(2)/u - t)/v))];

            ax = nexttile(tl, (iS-1)*nT + iT);
            histogram(ax, xData, 'BinEdges', edges, 'Normalization','pdf', ...
                'FaceColor',[0.7 0.7 0.7], 'EdgeColor','none');
            hold(ax,'on');
            yGrid = (0:length(pObs)-1)'*u;
            stairs(ax, yGrid, pObs/u, 'r', 'LineWidth', 1.5);
            hold(ax,'off');
            xlim(ax, xLim);
            if iS == 1
                title(ax, sprintf('t = %g min', times(iT)));
            end
            if iT == 1
                ylabel(ax, {sprintf('%s',Model.species{iS}),'density'});
            end
            if iS == nS
                xlabel(ax, 'intensity (a.u.)');
            end
        end
    end
    legend(nexttile(tl,1), {'data','model + PDO'}, 'Location','northeast');
end
end
