function figHandles = makeNucCytScatterPlots(ssaSoln,extendedMod,fignum,speciesIndMod, ...
    speciesIndDat,splitReps,crnadat,ctime,creplica,...
    transformSSAsoln)
arguments
    ssaSoln
    extendedMod
    fignum = 1;
    speciesIndMod = [1,2];
    speciesIndDat = [1,2];
    splitReps = false;
    crnadat = [3,4];
    ctime = 13;
    creplica = [];
    transformSSAsoln=[];
end

% Determine whether replica splitting is possible.
% Newer SSIT data sets store replica/meta-data columns separately in
% extendedMod.dataSet.savedData, with one row per row of extendedMod.dataSet.DATA.
% Older data sets may still store replica labels directly in DATA.
if splitReps
    if isempty(creplica)
        % Default to a saved column named 'replica' if available.
        if isfield(extendedMod.dataSet,'savedData') && ...
                istable(extendedMod.dataSet.savedData) && ...
                any(strcmp(extendedMod.dataSet.savedData.Properties.VariableNames,'replica'))
            creplica = 'replica';
        else
            splitReps = false;
        end
    elseif ischar(creplica) || isstring(creplica)
        creplica = char(creplica);

        if ~isfield(extendedMod.dataSet,'savedData') || ...
                ~istable(extendedMod.dataSet.savedData)
            error('makeNucCytScatterPlots:MissingSavedData', ...
                'splitReps=true, but extendedMod.dataSet.savedData is missing.');
        end

        if ~any(strcmp(extendedMod.dataSet.savedData.Properties.VariableNames,creplica))
            error('makeNucCytScatterPlots:MissingReplicaColumn', ...
                'Replica column "%s" was not found in extendedMod.dataSet.savedData.', ...
                creplica);
        end

    elseif isnumeric(creplica)
        % If a legacy numeric DATA column index was supplied but the new
        % savedData.replica column exists, prefer the new savedData path
        % whenever the numeric column is not present in DATA.
        if max(creplica) > size(extendedMod.dataSet.DATA,2) && ...
                isfield(extendedMod.dataSet,'savedData') && ...
                istable(extendedMod.dataSet.savedData) && ...
                any(strcmp(extendedMod.dataSet.savedData.Properties.VariableNames,'replica'))
            creplica = 'replica';
        elseif max(creplica) > size(extendedMod.dataSet.DATA,2)
            error('makeNucCytScatterPlots:ReplicaIndexOutOfBounds', ...
                ['Requested replica column index exceeds size(extendedMod.dataSet.DATA,2). ', ...
                 'Use a savedData column name such as ''replica'' instead.']);
        end
    end
end

% Transformation to results of SSA.  For example, could combine two species
% into one.
if ~isempty(transformSSAsoln)
    ssaSoln.trajs = tensorprod(transformSSAsoln, ssaSoln.trajs, 2, 1);
end


figHandles = figure;
    timeIndsDat = [1:length(extendedMod.dataSet.times)];

times2plot = extendedMod.dataSet.times(timeIndsDat);

Nrows = ceil(sqrt(length(timeIndsDat)));
Ncols = ceil(length(timeIndsDat)/Nrows);
for i = 1:length(timeIndsDat)
    subplot(Ncols,Nrows,i)

    % Find time in SSA data
    time = times2plot(i);
    [~,jSp] = min(abs(time-ssaSoln.T_array));
    
    % Add SSA to histogram plot
    M = squeeze(ssaSoln.trajs(speciesIndMod,jSp,:));
    % scatter(M(1,:),M(2,:),'ro');
    ksdensity(M','PlotFcn','contour');
    hold on

    % Add data to scatter plot
    timeMask = [extendedMod.dataSet.DATA{:,ctime}] == time;

    if splitReps
        if ischar(creplica) || isstring(creplica)
            % New preferred behavior:
            % replica/meta-data are stored in extendedMod.dataSet.savedData,
            % row-aligned with extendedMod.dataSet.DATA after filtering.
            replicaVals = extendedMod.dataSet.savedData.(char(creplica));
            replicaVals = replicaVals(timeMask);

            dMat = extendedMod.dataSet.DATA(timeMask,crnadat);
            repNames = unique(replicaVals);

            for j = 1:length(repNames)
                if iscell(replicaVals) || isstring(replicaVals) || iscategorical(replicaVals)
                    repMask = strcmp(string(replicaVals), string(repNames(j)));
                else
                    repMask = replicaVals == repNames(j);
                end

                dMatB = cell2mat(dMat(repMask,:));
                scatter(dMatB(:,1),dMatB(:,2),80,'.');
                hold on
            end

        else
            % Backward-compatible behavior:
            % replica labels were stored directly in extendedMod.dataSet.DATA.
            dMat = extendedMod.dataSet.DATA(timeMask,[crnadat,creplica]);
            repNames = unique(dMat(:,end));

            for j = 1:length(repNames)
                repMask = strcmp(string(dMat(:,end)), string(repNames(j)));
                dMatB = cell2mat(dMat(repMask,1:numel(crnadat)));
                scatter(dMatB(:,1),dMatB(:,2),80,'.');
                hold on
            end
        end
    else
        dMatB = cell2mat(extendedMod.dataSet.DATA(timeMask,crnadat));
        scatter(dMatB(:,1),dMatB(:,2),80,'k.');
        hold on
    end
    xyAll = cell2mat(extendedMod.dataSet.DATA(timeMask,crnadat));

    % dMat = extendedMod.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor;
    % 
    % xy = dMat.subs(dMat.subs(:,1)==timeIndsDat(i),1+speciesIndDat);   
    % v = dMat.values(dMat.subs(:,1)==timeIndsDat(i));
    % xyAll = zeros(sum(v),2);
    % for j = 1:length(v)
    %     xyAll(sum(v(1:j-1))+1:sum(v(1:j)),:) = repmat(xy(j,:),v(j),1);
    % end    
    % scatter(xyAll(:,1),xyAll(:,2),20,'k.');

    par0 = mean(M');
    cov12 = cov(M'); 
    ssit.parest.ellipse(par0,icdf('chi2',0.68,2)*cov12,'b','linewidth',3);  hold on;

    par0 = mean(xyAll);
    cov12 = cov(xyAll); 
    ssit.parest.ellipse(par0,icdf('chi2',0.68,2)*cov12,'c--','linewidth',3);  hold on;
    
    title(['t = ',num2str(time),'min']);
    xlabel('Nuc');ylabel('Cyt');
    set(gca,'xlim',[0,200],'ylim',[0,300])
    
end
end
