function [T_all, figHandles] = fanosFromSSATrajs(modelSpecs, expDataSet, speciesIdxModel, speciesLabel, speciesIdxData, opts)
% fanoFromSSATrajs
%
% Computes means, variances, and Fano factors from an arbitrary number of
% SSIT SSA trajectory solutions, compares them to an experimental SSIT
% dataSet, and generates four plots:
%
%   1. Variance vs. mean
%   2. Distribution of selected species at one specified time point
%   3. Mean vs. time
%   4. Variance vs. time
%
% Expected model solution format:
%   Model.Solutions.trajs   = nSpecies x nTimes x nTrajs
%   Model.Solutions.T_array = 1 x nTimes
%
% Inputs:
%   modelSpecs  - struct array with fields:
%                   .Model
%                   .Label
%
%                 Optional fields:
%                   .LineSpec
%
%   expDataSet  - SSIT dataSet struct with fields:
%                   .times
%                   .mean
%                   .var
%                   .nCells
%
%   speciesIdxModel  - index of species in SSA trajectories
%   speciesLabel - display name, e.g. "mRNA_{cyt}"
%
% Options:
%   opts.TimePoint         - time point for distribution plot, default 75
%   opts.ExcludeTimesBelow - exclude negative initialization times, default 0
%   opts.MakePlots         - true/false, default true
%   opts.SaveFigures       - true/false, default false
%   opts.OutDir            - output directory, default "savedFigures"
%   opts.FilePrefix        - prefix for saved figures
%   opts.HistBinMethod     - histogram bin method, default "integers"
%   opts.UseNearestTime    - use nearest time if exact time absent, default true
%
% Outputs:
%   T_all      - long-format table with model and data moments
%   figHandles - struct of figure handles

arguments
    modelSpecs struct
    expDataSet struct
    speciesIdxModel double 
    speciesLabel string
    speciesIdxData double 
    opts.TimePoint (1,1) double = 75
    opts.ExcludeTimesBelow (1,1) double = 0
    opts.MakePlots logical = true
    opts.SaveFigures logical = false
    opts.OutDir string = "savedFigures"
    opts.FilePrefix string = "FanoComparison"
    opts.HistBinMethod string = "integers"
    opts.UseNearestTime logical = true
    opts.showErrorBars = false
end

figHandles = struct();

%% ------------------------------------------------------------------------
%  Validate modelSpecs
% -------------------------------------------------------------------------

if ~isfield(modelSpecs, "Model")
    error('modelSpecs must contain a field named Model.')
end

if ~isfield(modelSpecs, "Label")
    error('modelSpecs must contain a field named Label.')
end

nModels = numel(modelSpecs);

%% ------------------------------------------------------------------------
%  Compute model moments
% -------------------------------------------------------------------------

T_model_all = table();

for iModel = 1:nModels
    M = modelSpecs(iModel).Model;
    label = string(modelSpecs(iModel).Label);

    if ~isfield(M.Solutions, "trajs")
        error('Model "%s" does not contain Solutions.trajs.', label)
    end

    if ~isfield(M.Solutions, "T_array")
        error('Model "%s" does not contain Solutions.T_array.', label)
    end

    trajs = M.Solutions.trajs;

    % Perform grouping of bound and unbound RNA for total Cyt count.
    if ~isempty(modelSpecs(iModel).transformSSAsoln)
        trajs = tensorprod(modelSpecs(iModel).transformSSAsoln, trajs, 2, 1);
    end

    times = M.Solutions.T_array(:);

    if speciesIdxModel > size(trajs, 1)
        error('speciesIdxModel=%d exceeds number of species in model "%s".', ...
            speciesIdxModel, label)
    end

    nTimes = numel(times);
    nTrajs = size(trajs, 3);

    mu = zeros(nTimes, 1);
    sig2 = zeros(nTimes, 1);
    fano = zeros(nTimes, 1);

    for it = 1:nTimes
        vals = squeeze(trajs(speciesIdxModel, it, :));
        vals = vals(:);

        mu(it) = mean(vals);
        sig2(it) = var(vals, 0);  
        fano(it) = sig2(it) ./ max(mu(it), eps);
    end

    T_tmp = table;
    T_tmp.Source = repmat(label, nTimes, 1);
    T_tmp.Type = repmat("Model", nTimes, 1);
    T_tmp.Species = repmat(speciesLabel, nTimes, 1);
    T_tmp.Time_min = times;
    T_tmp.N = repmat(nTrajs, nTimes, 1);
    T_tmp.Mean = mu;
    T_tmp.Var = sig2;
    T_tmp.Fano = fano;

    T_tmp.MeanStd = NaN(nTimes, 1);
    T_tmp.VarStd = NaN(nTimes, 1);
    T_tmp.FanoStd = NaN(nTimes, 1);

    T_model_all = [T_model_all; T_tmp]; 
end

%% ------------------------------------------------------------------------
%  Experimental moments table
% -------------------------------------------------------------------------

%% Extract replica data
dataFileDusp1 = 'Data/DUSP1_SSITcellresults_Final_Sep18.csv';
allData = readtable(dataFileDusp1);
allData = allData(allData.cyto_area_px>=12593&allData.cyto_area_px<=17685&isnan(allData.time_tpl),:);
allData = allData((allData.dex_conc == 100)|(allData.time == 0),:);
%
means_v_t = cell(1,length(times));
for it = 1:length(times)
    time = times(it);
    allDataT = allData(allData.time == time,:);
    
    reps = unique(allDataT.replica);
    for j = 1:length(reps)
        means_v_t{it}.nuc(j) = mean(allDataT.num_nuc_noTS(strcmp(allDataT.replica,reps{j})));
        means_v_t{it}.cyt(j) = mean(allDataT.num_cyto_spots(strcmp(allDataT.replica,reps{j})));
        vars_v_t{it}.nuc(j) = var(allDataT.num_nuc_noTS(strcmp(allDataT.replica,reps{j})));
        vars_v_t{it}.cyt(j) = var(allDataT.num_cyto_spots(strcmp(allDataT.replica,reps{j})));
        fanos_v_t{it}.nuc(j) = vars_v_t{it}.nuc(j)/means_v_t{it}.nuc(j);
        fanos_v_t{it}.cyt(j) = vars_v_t{it}.cyt(j)/means_v_t{it}.cyt(j);
    end
    if ~isempty(means_v_t{it})
        mean_means_v_t(:,it) = [mean(means_v_t{it}.nuc);mean(means_v_t{it}.cyt)];
        std_means_v_t(:,it) = [std(means_v_t{it}.nuc);std(means_v_t{it}.cyt)];
        mean_vars_v_t(:,it) = [mean(vars_v_t{it}.nuc);mean(vars_v_t{it}.cyt)];
        std_vars_v_t(:,it) = [std(vars_v_t{it}.nuc);std(vars_v_t{it}.cyt)];
        mean_fanos_v_t(:,it) = [mean(fanos_v_t{it}.nuc);mean(fanos_v_t{it}.cyt)];
        std_fanos_v_t(:,it) = [std(fanos_v_t{it}.nuc);std(fanos_v_t{it}.cyt)];
    end
end

%% Convert replicate-based experimental summaries into table

% speciesIdxData: 1 = nuc, 2 = cyt, according to mean_means_v_t rows.
if speciesIdxData == 1
    repSpeciesField = "nuc";
elseif speciesIdxData == 2
    repSpeciesField = "cyt";
else
    error('speciesIdxData must be 1 for nuc or 2 for cyt for replicate summaries.')
end

validRepTimes = false(length(times),1);

for it = 1:length(times)
    validRepTimes(it) = exist('means_v_t','var') && ...
        length(means_v_t) >= it && ...
        ~isempty(means_v_t{it}) && ...
        isfield(means_v_t{it}, repSpeciesField);
end

repTimes = times(validRepTimes);

repMeanMean = mean_means_v_t(speciesIdxData, validRepTimes)';
repStdMean  = std_means_v_t(speciesIdxData, validRepTimes)';

repMeanVar = mean_vars_v_t(speciesIdxData, validRepTimes)';
repStdVar  = std_vars_v_t(speciesIdxData, validRepTimes)';

repMeanFano = mean_fanos_v_t(speciesIdxData, validRepTimes)';
repStdFano  = std_fanos_v_t(speciesIdxData, validRepTimes)';

T_data_rep = table;
T_data_rep.Source = repmat("Experimental data", numel(repTimes), 1);
T_data_rep.Type = repmat("Data", numel(repTimes), 1);
T_data_rep.Species = repmat(speciesLabel, numel(repTimes), 1);
T_data_rep.Time_min = repTimes(:);
T_data_rep.Mean = repMeanMean(:);
T_data_rep.MeanStd = repStdMean(:);
T_data_rep.Var = repMeanVar(:);
T_data_rep.VarStd = repStdVar(:);
T_data_rep.Fano = repMeanFano(:);
T_data_rep.FanoStd = repStdFano(:);

%%
dataTimes = expDataSet.times(:);

dataMean = expDataSet.mean(:,speciesIdxData);
dataVar  = expDataSet.var(:,speciesIdxData);

% Handle nCells
if isfield(expDataSet, "nCells")
    dataNRaw = expDataSet.nCells;

    if isvector(dataNRaw)
        dataN = dataNRaw(:);
    else
        if size(dataNRaw, 1) == numel(dataTimes)
            dataN = dataNRaw(:, 1);
        elseif size(dataNRaw, 2) == numel(dataTimes)
            dataN = dataNRaw(1, :)';
        else
            dataN = NaN(numel(dataTimes), 1);
        end
    end
else
    dataN = NaN(numel(dataTimes), 1);
end

% Defensive check
if numel(dataMean) ~= numel(dataTimes)
    error('Length mismatch: expDataSet.mean has %d values, but expDataSet.times has %d values.', ...
        numel(dataMean), numel(dataTimes))
end

if numel(dataVar) ~= numel(dataTimes)
    error('Length mismatch: expDataSet.var has %d values, but expDataSet.times has %d values.', ...
        numel(dataVar), numel(dataTimes))
end

if numel(dataN) ~= numel(dataTimes)
    dataN = NaN(numel(dataTimes), 1);
end

T_data = table;
T_data.Source = repmat("Experimental data", numel(dataTimes), 1);
T_data.Type = repmat("Data", numel(dataTimes), 1);
T_data.Species = repmat(speciesLabel, numel(dataTimes), 1);
T_data.Time_min = dataTimes;
T_data.N = dataN;
T_data.Mean = dataMean(:);
T_data.Var = dataVar(:);
T_data.Fano = dataVar(:) ./ max(dataMean(:), eps);

% Attach replicate-based uncertainty where available.
T_data.MeanStd = NaN(height(T_data),1);
T_data.VarStd = NaN(height(T_data),1);
T_data.FanoStd = NaN(height(T_data),1);

if exist('T_data_rep','var') && ~isempty(T_data_rep)
    [tfRep, idxRep] = ismember(T_data.Time_min, T_data_rep.Time_min);

    T_data.MeanStd(tfRep) = T_data_rep.MeanStd(idxRep(tfRep));
    T_data.VarStd(tfRep) = T_data_rep.VarStd(idxRep(tfRep));
    T_data.FanoStd(tfRep) = T_data_rep.FanoStd(idxRep(tfRep));
end

%% ------------------------------------------------------------------------
%  Combine and filter
% -------------------------------------------------------------------------

T_all = [T_model_all; T_data];

if ~isnan(opts.ExcludeTimesBelow)
    T_all = T_all(T_all.Time_min >= opts.ExcludeTimesBelow, :);
end

%% ------------------------------------------------------------------------
%  Plotting
% -------------------------------------------------------------------------

if opts.MakePlots

    if opts.SaveFigures && ~exist(opts.OutDir, "dir")
        mkdir(opts.OutDir);
    end

    %% 1. Variance vs mean
    figHandles.varVsMean = figure;
    hold on

    for iModel = 1:nModels
        label = string(modelSpecs(iModel).Label);
        idx = T_all.Source == label;

        plot(T_all.Mean(idx), T_all.Var(idx), ...
            '-o', 'LineWidth', 2, 'DisplayName', label);
    end

    % idxData = T_all.Type == "Data";
    % plot(T_all.Mean(idxData), T_all.Var(idxData), ...
    %     'ko', 'MarkerFaceColor', 'k', ...
    %     'LineStyle', '--', ...
    %     'DisplayName', 'Experimental data');

    idxData = T_all.Type == "Data";

    xData = T_all.Mean(idxData);
    yData = T_all.Var(idxData);
    xErr  = T_all.MeanStd(idxData);
    yErr  = T_all.VarStd(idxData);
    
    if opts.showErrorBars
        plotErrorbarXY(xData, yData, xErr, yErr, ...
            'ko', 'MarkerFaceColor', 'k', ...
            'LineStyle', '--', ...
            'LineWidth', 1.5, ...
            'DisplayName', 'Experimental data');
    else
        plot(xData, yData, ...
            'ko', 'MarkerFaceColor', 'k', ...
            'LineStyle', '--', ...
            'LineWidth', 1.5, ...
            'DisplayName', 'Experimental data');
    end

    xlabel(['Mean ', char(speciesLabel)], 'Interpreter', 'tex', FontSize=24)
    ylabel(['Variance ', char(speciesLabel)], 'Interpreter', 'tex', FontSize=24)
    title(['Variance vs. mean for ', char(speciesLabel)], 'Interpreter', 'tex', FontSize=24)
    legend('Location', 'northwest', FontSize=16)
    set(gca, 'FontSize', 20)
    box on
    grid on

    if opts.SaveFigures
        savefig(figHandles.varVsMean, fullfile(opts.OutDir, opts.FilePrefix + "_VarVsMean.fig"));
        saveas(figHandles.varVsMean, fullfile(opts.OutDir, opts.FilePrefix + "_VarVsMean.png"));
    end

    %% 2. Distribution at specified time
    figHandles.distribution = figure;
    hold on

    for iModel = 1:nModels
        M = modelSpecs(iModel).Model;
        label = string(modelSpecs(iModel).Label);

        times = M.Solutions.T_array(:);
        it = findTimeIndex(times, opts.TimePoint, opts.UseNearestTime);

        % Group lumped species into one.
        trajs = M.Solutions.trajs;
        if ~isempty(modelSpecs(iModel).transformSSAsoln)
            trajs = tensorprod(modelSpecs(iModel).transformSSAsoln, trajs, 2, 1);
        end

        vals = squeeze(trajs(speciesIdxModel, it, :));
        vals = vals(:);

        histogram(vals, ...
            'Normalization', 'pdf', ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2, ...
            'BinMethod', opts.HistBinMethod, ...
            'DisplayName', label);
    end

    % Try to overlay experimental single-cell data at this time, if raw data
    % can be extracted from expDataSet.DATA.
    [dataVals, gotDataVals] = tryExtractExperimentalValuesAtTime(expDataSet, opts.TimePoint);

    if gotDataVals
        histogram(dataVals(:,speciesIdxData), ...
            'Normalization', 'pdf', ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2, ...
            'DisplayName', 'Experimental data',...
            'EdgeColor','k');
    else
        warning(['Could not extract raw experimental values for distribution plot. ', ...
                 'Only model distributions were plotted.'])
    end

    xlabel(char(speciesLabel), 'Interpreter', 'tex')
    ylabel('Probability')
    title(sprintf('Distribution of %s at t = %.3g min', char(speciesLabel), opts.TimePoint), ...
        'Interpreter', 'tex')
    legend('Location', 'northeast', FontSize=16)
    set(gca, 'FontSize', 20)
    box on
    grid on

    if opts.SaveFigures
        savefig(figHandles.distribution, fullfile(opts.OutDir, opts.FilePrefix + "_Distribution_t" + num2str(opts.TimePoint) + ".fig"));
        saveas(figHandles.distribution, fullfile(opts.OutDir, opts.FilePrefix + "_Distribution_t" + num2str(opts.TimePoint) + ".png"));
    end

    %% 3. Mean vs time
    figHandles.meanVsTime = figure;
    hold on

    for iModel = 1:nModels
        label = string(modelSpecs(iModel).Label);
        idx = T_all.Source == label;

        plot(T_all.Time_min(idx), T_all.Mean(idx), ...
            '-o', 'LineWidth', 2, 'DisplayName', label);
    end

    if ~opts.showErrorBars
        plot(T_data.Time_min, T_data.Mean, ...
            'ko', 'MarkerFaceColor', 'k', ...
            'LineStyle', '--', ...
            'DisplayName', 'Experimental data');

    else
        errorbar(T_data.Time_min, T_data.Mean, T_data.MeanStd, ...
            'ko', ...
            'MarkerFaceColor', 'k', ...
            'LineStyle', '--', ...
            'LineWidth', 1.5, ...
            'CapSize', 8, ...
            'DisplayName', 'Experimental data');
    end

    xlabel('Time (min)', FontSize=24)
    ylabel(['Mean ', char(speciesLabel)], 'Interpreter', 'tex', FontSize=24)
    title(['Mean vs. time for ', char(speciesLabel)], 'Interpreter', 'tex', FontSize=24)
    legend('Location', 'southeast', FontSize=16)
    set(gca, 'FontSize', 20)
    box on
    grid on

    if opts.SaveFigures
        savefig(figHandles.meanVsTime, fullfile(opts.OutDir, opts.FilePrefix + "_MeanVsTime.fig"));
        saveas(figHandles.meanVsTime, fullfile(opts.OutDir, opts.FilePrefix + "_MeanVsTime.png"));
    end

    %% 4. Variance vs time
    figHandles.varVsTime = figure;
    hold on

    for iModel = 1:nModels
        label = string(modelSpecs(iModel).Label);
        idx = T_all.Source == label;

        plot(T_all.Time_min(idx), T_all.Var(idx), ...
            '-o', 'LineWidth', 2, 'DisplayName', label);
    end

    if ~opts.showErrorBars
        plot(T_data.Time_min, T_data.Var, ...
            'ko', 'MarkerFaceColor', 'k', ...
            'LineStyle', '--', ...
            'DisplayName', 'Experimental data');
    else
        errorbar(T_data.Time_min, T_data.Var, T_data.VarStd, ...
            'ko', ...
            'MarkerFaceColor', 'k', ...
            'LineStyle', '--', ...
            'LineWidth', 1.5, ...
            'CapSize', 8, ...
            'DisplayName', 'Experimental data');
    end

    xlabel('Time (min)', FontSize=24)
    ylabel(['Variance ', char(speciesLabel)], 'Interpreter', 'tex', FontSize=24)
    title(['Variance vs. time for ', char(speciesLabel)], 'Interpreter', 'tex', FontSize=24)
    legend('Location', 'southeast', FontSize=16)
    set(gca, 'FontSize', 20)
    box on
    grid on

    if opts.SaveFigures
        savefig(figHandles.varVsTime, fullfile(opts.OutDir, opts.FilePrefix + "_VarVsTime.fig"));
        saveas(figHandles.varVsTime, fullfile(opts.OutDir, opts.FilePrefix + "_VarVsTime.png"));
    end
end

end

%% ========================================================================
%  Local helper: 2D error bars
% ========================================================================

function h = plotErrorbarXY(x, y, xErr, yErr, varargin)
% plotErrorbarXY
% Draws symmetric horizontal and vertical error bars.
%
% Inputs:
%   x, y   - point coordinates
%   xErr   - symmetric x uncertainty
%   yErr   - symmetric y uncertainty
%
% Usage:
%   plotErrorbarXY(x, y, xErr, yErr, 'ko', 'MarkerFaceColor', 'k')

x = x(:);
y = y(:);
xErr = xErr(:);
yErr = yErr(:);

holdState = ishold;
hold on

% Remove entries with missing central values.
valid = ~isnan(x) & ~isnan(y);
x = x(valid);
y = y(valid);
xErr = xErr(valid);
yErr = yErr(valid);

% Replace missing errors with zero.
xErr(isnan(xErr)) = 0;
yErr(isnan(yErr)) = 0;

% Main points.
h = plot(x, y, varargin{:});

% Error bars should not clutter legend.
for i = 1:numel(x)
    plot([x(i)-xErr(i), x(i)+xErr(i)], [y(i), y(i)], ...
        'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');

    plot([x(i), x(i)], [y(i)-yErr(i), y(i)+yErr(i)], ...
        'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
end

if ~holdState
    hold off
end

end


%% ========================================================================
%  Local helper: find requested time index
% ========================================================================

function it = findTimeIndex(times, requestedTime, useNearest)

times = times(:);

idx = find(abs(times - requestedTime) < 1e-9, 1);

if ~isempty(idx)
    it = idx;
    return
end

if useNearest
    [~, it] = min(abs(times - requestedTime));
    warning('Requested time %.3g not found exactly. Using nearest time %.3g.', ...
        requestedTime, times(it));
else
    error('Requested time %.3g not found in T_array.', requestedTime)
end

end


%% ========================================================================
%  Local helper: try to extract experimental single-cell values
% ========================================================================

function [vals, success] = tryExtractExperimentalValuesAtTime(ds, requestedTime)

% success = false;
% vals = [];
% 
% if ~isfield(ds, "DATA") || isempty(ds.DATA)
%     return
% end
% DATA = ds.DATA;

t = [ds.DATA{:,1}]';
vals = cell2mat(ds.DATA(t==requestedTime,2:end));
success = true;

% % Common case: DATA is an nCells x 2 cell array where one column is time
% % and the other column is the observed value.
% if iscell(DATA) && size(DATA, 2) >= 2
% 
%     try
%         col1 = DATA(:,2);
%         col2 = DATA(:,3);
% 
%         t = [DATA{:,1}]';
%         x1 = cellToNumericVector(col1);
%         x2 = cellToNumericVector(col2);
% 
%         % Heuristic: whichever column contains requestedTime and resembles
%         % the dataSet.times vector is treated as time.
%         if any(abs(t - requestedTime) < 1e-9)
%             timeVals = x1;
%             obsVals = x2;
%         elseif any(abs(x2 - requestedTime) < 1e-9)
%             timeVals = x2;
%             obsVals = x1;
%         else
%             return
%         end
% 
%         idx = abs(timeVals - requestedTime) < 1e-9;
%         vals = obsVals(idx);
%         vals = vals(~isnan(vals));
% 
%         if ~isempty(vals)
%             success = true;
%         end
% 
%     catch
%         success = false;
%         vals = [];
%     end
% end

end


%% ========================================================================
%  Local helper: cell column to numeric vector
% ========================================================================

function x = cellToNumericVector(c)

x = NaN(numel(c), 1);

for i = 1:numel(c)
    ci = c{i};

    if isnumeric(ci) && isscalar(ci)
        x(i) = ci;
    elseif islogical(ci) && isscalar(ci)
        x(i) = double(ci);
    elseif ischar(ci) || isstring(ci)
        tmp = str2double(ci);
        if ~isnan(tmp)
            x(i) = tmp;
        end
    end
end

end