function [KS,figHandles] = makeCytDistPlots(ssaSoln,model,fignum,speciesIndMod, ...
    speciesIndDat,binEdges,splitReps,plotType,crnadat,ctime,creplica,quartiles, ...
    transformSSAsoln)
arguments
    ssaSoln
    model
    fignum = 1;
    speciesIndMod = 1;
    speciesIndDat = 1;
    binEdges = [0:10:300];
    splitReps = false
    plotType = 'pdf'
    crnadat = [3,4];
    ctime = [13];
    creplica = [];
    quartiles = true;
    transformSSAsoln = [];
end

% if isempty(creplica)
%     splitReps = false;
% end

% Determine whether replica splitting is possible.
% In newer SSIT data sets, replica/meta-data columns are stored separately
% in model.dataSet.savedData, with column names in
% model.dataSet.savedColumns:
if splitReps
    if isempty(creplica)
        % Default to a saved column named 'replica' if available.
        if isfield(model.dataSet,'savedData') && ...
                istable(model.dataSet.savedData) && ...
                any(strcmp(model.dataSet.savedData.Properties.VariableNames,'replica'))
            creplica = 'replica';
        else
            splitReps = false;
        end
    elseif ischar(creplica) || isstring(creplica)
        creplica = char(creplica);

        if ~isfield(model.dataSet,'savedData') || ~istable(model.dataSet.savedData)
            error('makeCytDistPlots:MissingSavedData', ...
                'splitReps=true, but model.dataSet.savedData is missing.');
        end

        if ~any(strcmp(model.dataSet.savedData.Properties.VariableNames,creplica))
            error('makeCytDistPlots:MissingReplicaColumn', ...
                'Replica column "%s" was not found in model.dataSet.savedData.', creplica);
        end
    end
end

% Transformation to results of SSA.  For example, could combine two species
% into one.
if ~isempty(transformSSAsoln)
    ssaSoln.trajs = tensorprod(transformSSAsoln, ssaSoln.trajs, 2, 1);
end

figHandles{1} = figure;
timeIndsDat = [1:length(model.dataSet.times)];

times2plot = model.dataSet.times;

Nrows = ceil(sqrt(length(timeIndsDat)));
Ncols = ceil(length(timeIndsDat)/Nrows);
KS = zeros(1,length(timeIndsDat));
for i = 1:length(timeIndsDat)
    subplot(Ncols,Nrows,i); hold off
    time = times2plot(i);
    
    % Make plots for DATA
    % if splitReps
    %     dMat = model.dataSet.DATA([model.dataSet.DATA{:,ctime}]==time,[crnadat,creplica]);
    %     repNames = unique(dMat(:,3));
    %     for j = 1:length(repNames)
    %         dMatB = cell2mat(dMat(strcmp(dMat(:,3),repNames{j}),1:2));
    %         h = histogram(dMatB(:,speciesIndDat),binEdges,'Normalization',plotType,'DisplayStyle','stairs','LineWidth',2);
    %         hold on
    %     end
    % else
    %     dMatB = cell2mat(model.dataSet.DATA([model.dataSet.DATA{:,ctime}]==time,crnadat));
    %     histogram(dMatB(:,speciesIndDat),binEdges,'Normalization',plotType,'DisplayStyle','stairs','LineWidth',3,'EdgeColor',[.40,.40,.80]);
    %     hold on
    % end
    % dMatB = cell2mat(model.dataSet.DATA([model.dataSet.DATA{:,ctime}]==time,crnadat));

    % Make plots for DATA
    timeMask = [model.dataSet.DATA{:,ctime}] == time;
    
    if splitReps
    
        % New preferred behavior:
        % replica/meta-data are stored in model.dataSet.savedData, row-aligned
        % with model.dataSet.DATA after filtering:
        if ischar(creplica) || isstring(creplica)
            replicaVals = model.dataSet.savedData.(char(creplica));
            replicaVals = replicaVals(timeMask);
    
            dMat = model.dataSet.DATA(timeMask,crnadat);
    
            repNames = unique(replicaVals);
    
            for j = 1:length(repNames)
                if iscell(replicaVals) || isstring(replicaVals) || iscategorical(replicaVals)
                    repMask = strcmp(string(replicaVals), string(repNames(j)));
                else
                    repMask = replicaVals == repNames(j);
                end
    
                dMatB_rep = cell2mat(dMat(repMask,:));
    
                histogram(dMatB_rep(:,speciesIndDat),binEdges, ...
                    'Normalization',plotType, ...
                    'DisplayStyle','stairs', ...
                    'LineWidth',2);
                hold on
            end
    
        else
            % Backward-compatible behavior:
            % replica was stored directly in model.dataSet.DATA:
            dMat = model.dataSet.DATA(timeMask,[crnadat,creplica]);
            repNames = unique(dMat(:,end));
    
            for j = 1:length(repNames)
                repMask = strcmp(string(dMat(:,end)), string(repNames{j}));
                dMatB_rep = cell2mat(dMat(repMask,1:numel(crnadat)));
    
                histogram(dMatB_rep(:,speciesIndDat),binEdges, ...
                    'Normalization',plotType, ...
                    'DisplayStyle','stairs', ...
                    'LineWidth',2);
                hold on
            end
        end
    
    else
        dMatB = cell2mat(model.dataSet.DATA(timeMask,crnadat));
    
        histogram(dMatB(:,speciesIndDat),binEdges, ...
            'Normalization',plotType, ...
            'DisplayStyle','stairs', ...
            'LineWidth',3, ...
            'EdgeColor',[.40,.40,.80]);
        hold on
    end
    
    % Always keep the full filtered data for KS/statistics below.
    dMatB = cell2mat(model.dataSet.DATA(timeMask,crnadat));

    % Add SSA to histogram plot
    % Find time in SSA data
    [~,jSp] = min(abs(time-ssaSoln.T_array));
    M = squeeze(ssaSoln.trajs(speciesIndMod,jSp,:));
    histogram(M,binEdges,'Normalization',plotType,'DisplayStyle','stairs','LineWidth',3,'EdgeColor',[0,0,0]);

    [~,~,KS(i)] = kstest2(dMatB(:,speciesIndDat),M');
    % % Add data to histogram plot
    % dMat = double(extendedMod.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor(timeIndsDat(i),:,:));
    % N = sum(dMat,"all");
    % if speciesIndDat==1
    %     PD = [0;sum(dMat,2)/N];
    % else
    %     PD = [0;sum(dMat,1)'/N];
    % end
    % binEdges = round(H.BinEdges);
    % PD(binEdges(end)+1) = 0;
    % 
    % nBins = length(binEdges)-1;
    % PDbinned = zeros(nBins+1,1);
    % binwidth = binEdges(2)-binEdges(1);
    % for j = 1:nBins
    %     PDbinned(j) = sum(PD(binEdges(j)+1:binEdges(j+1)));
    % end
    % 
    % PDbinned(end) = 1-sum(PDbinned);
    % stairs(binEdges,PDbinned/binwidth,'linewidth',2)
    title(['t = ',num2str(time),'min']);

    % set(gca,'FontSize',15,'ylim',[0,0.03])

    
    modMean(i) = mean(M);
    modSTD(i) = std(M);
    mSorted = sort(M); l = length(M);
    modQuartiles(:,i) = mSorted([floor(l/4),floor(l/2),ceil(l*3/4)]);
   
    datMean(i) = mean(dMatB(:,speciesIndDat));
    datSTD(i) = std(dMatB(:,speciesIndDat));
    dSorted = sort(dMatB(:,speciesIndDat)); l = length(dSorted);
    datQuartiles(:,i) = dSorted([floor(l/4),floor(l/2),ceil(l*3/4)]);
    
end

figHandles{2} = figure;

vx = [times2plot,times2plot(end:-1:1)];
if quartiles
    vy = [modQuartiles(1,:),modQuartiles(3,end:-1:1)];
    fi = fill(vx, vy, 'b'); hold on;
    plot(times2plot,modQuartiles(2,:),'b','LineWidth',2);
    eb = errorbar(times2plot,datQuartiles(2,:),datQuartiles(2,:)-datQuartiles(1,:),datQuartiles(3,:)-datQuartiles(2,:));
else
    vy = [modMean-modSTD,modMean(end:-1:1)+modSTD(end:-1:1)];
    fi = fill(vx, vy, 'b'); hold on;
    plot(times2plot,modMean,'b','LineWidth',2);
    eb = errorbar(times2plot,datMean,datSTD);
end

fi.FaceColor = [0.9,0.9,1];

eb.LineWidth = 2;
eb.Color = 'b';
eb.Marker = 'o';
eb.MarkerSize = 8;
eb.MarkerFaceColor = 'k';

end
