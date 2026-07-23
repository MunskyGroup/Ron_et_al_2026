
function [TSLikelihood,modelResults,dataResults] = computeTSlikelihood(pars, ...
    ModelTS,allData,outlierThresh,reportErrors)
arguments
    pars
    ModelTS
    allData
    outlierThresh = 50;  
    reportErrors = false
    % The outlierThresh input sets an outlier threshold for the nascent RNA
    % at the TS.  Data with observed counts more than this value will be
    % trincated to this value.  Model predictions for counts higher than
    % this value will also be integrated into a single bin "more than
    % outlierThresh".
end
    
ModelTS.parameters(ModelTS.fittingOptions.modelVarsToFit,2) = num2cell(pars);
tau_elong = ModelTS.parameters{19,2};
ModelTS.parameters{19,2} = 0; % tau elong is not used as a delay in the TS model.
dex = ModelTS.parameters{13,2};

% Set the min mumber of nascent RNA to be considerred a TS. this should
% match the value set in the image processing.
TSthresh = 4;

% The rest of the parameters go into the model
% x = x(1:end-1);
% ModelTS.parameters(ModelTS.fittingOptions.modelVarsToFit,2) = num2cell(x);
% ModelTS.parameters(13,:) = {'Dex0',dex};

% set upper bound on FSP solution
nSp = length(ModelTS.species) - length(ModelTS.hybridOptions.upstreamODEs);
ModelTS.fspOptions.bounds(nSp+3) = outlierThresh*1.5;

% turn off FSP expansion routine
ModelTS.fspOptions.fspTol = inf;

% Create past model only used to get distribution of states for when
% elongation would begin.
ModelTS.fspOptions.stateSpace = [];
ModelTSPast = ModelTS;
ModelTSPast.parameters(3,:) = {'kr_on_0',0};  % Turn off transcription.
ModelTSPast.parameters(23,:) = {'kr_off',0};  % Turn off transcription.
ModelTSPast.parameters(26,:) = {'kr_on_1',0};  % Turn off transcription.

% For TS calculations, we consider only one alele at a time.
ModelTSPast.initialCondition(1) = 1;
ModelTSPast.initialCondition(5) = 0;
ModelTSPast.fspOptions.initApproxSS = false;

% Shift solution time points back by tau_elong.
ModelTSPast.tSpan = ModelTS.tSpan - tau_elong;
ModelTSPast.tSpan = [-200-tau_elong,ModelTSPast.tSpan];
ModelTSPast.initialTime = -200-tau_elong;

% Solve for distributions at t - tau_elong.
% if nargout==1
%     ModelTSPast.tSpan = ModelTSPast.tSpan(1:2);
% end
TSPastSoln = ModelTSPast.solve(returnType='soln');
NT = length(TSPastSoln.fsp)-1;

% For the present model we start one elongation period in the past and then 
% integrate forward in time.
ModelTSPresent = ModelTS;
% During transcription, there is no degradation/transport.
ModelTSPresent.parameters(4,:) = {'knuc2cyt',0};

% We will NOT use steady state initial conditions -- rather, we use the
% computed distributions from previous time step.
ModelTSPresent.fspOptions.initApproxSS = false;
% We need to solve once for each time point, where the initial condition is
% defined by the past model.

fsp = cell(1,NT);
for iT = 1:NT
    tmpModelTSPresent = ModelTSPresent.setICfromFspVector(TSPastSoln.stateSpace,TSPastSoln.fsp{iT+1});
    tmpModelTSPresent.tSpan = ModelTSPast.tSpan(iT+1)+[0,tau_elong/2,tau_elong];
    tmpModelTSPresent.initialTime = tmpModelTSPresent.tSpan(1);
    tsPresentSoln = tmpModelTSPresent.solve(returnType='soln');
    fsp(iT) = tsPresentSoln.fsp(end); 
end

% Summarize the model results for plotting and for likelihood function
% calculations.
modelResults.meanNascentMod = zeros(NT,1);
modelResults.fracTSMod = zeros(NT,1);
for iT = NT:-1:1
    if fsp{iT}.p.dim==3
        P=double(fsp{iT}.p.sumOver([1,2]).data);
    else
        P=double(fsp{iT}.p.sumOver([1,2,4]).data);
    end
  
    P = max(0,P);
    % if sum(P)>1
    P = P/sum(P);
    % end
    % P(1:TSthresh) = [sum(P(1:TSthresh));zeros(TSthresh-1,1)]; 
    % modelResults.Psaved{iT} = P;
    modelResults.Psaved{iT} = max(1e-10,P);
end

% Summarize the data for plotting and for likelihood function
% calculations.
for iT = NT:-1:1
    time = ModelTS.tSpan(iT);
    inds = (allData.dex_conc==dex&allData.time==time)|(allData.time==time&time==0);
    redData = allData(inds,:);

    % remove outliers with counts above threshold
    redData = redData(redData.largest_ts<=outlierThresh,:);

    TS_counts0 = redData.largest_ts;
    TS_counts1 = redData.second_largest_ts;    
    TS_counts0(isnan(TS_counts0)) = 0;
    TS_counts1(isnan(TS_counts1)) = 0;
    dataResults.TS_counts{iT} = TS_counts0+TS_counts1; 
    dataResults.meanNascentDat(iT) = mean(dataResults.TS_counts{iT});
    isTS = dataResults.TS_counts{iT}>0;
    dataResults.meanNascentTSDat(iT) = mean(dataResults.TS_counts{iT}(isTS));
    dataResults.fracTSDat(iT) = sum(dataResults.TS_counts{iT}>=TSthresh)/length(dataResults.TS_counts{iT});
    dataResults.N(iT) = sum(isTS);
    dataResults.meanNascentTSDatstd(iT) = std(dataResults.TS_counts{iT}(isTS))./sqrt(dataResults.N(iT));

    if ~isempty(isTS)&dataResults.N(iT)==0
        dataResults.meanNascentDat(iT) = 0;
        dataResults.meanNascentTSDat(iT) = 0;
        dataResults.fracTSDat(iT) = 0;
        dataResults.fracTSDatstd(iT) = 0;
        dataResults.meanNascentTSDatstd(iT) = 0;
        dataResults.meanNascentDatstd(iT) = 0;
        dataResults.N(iT)=1;
    end

    if reportErrors
        reps = unique(redData.replica);
        frac = [];
        mn = [];
        mnTS = [];
        for j = 1:length(reps)
            AllDatRedReps = redData(strcmp(redData.replica,reps{j}),:);
            frac(j) = sum(AllDatRedReps.largest_ts>0)/size(AllDatRedReps,1);
            mn(j) = (sum(AllDatRedReps.largest_ts(AllDatRedReps.largest_ts>0))+...
                sum(AllDatRedReps.second_largest_ts(AllDatRedReps.second_largest_ts>0)))/size(AllDatRedReps,1);
            isTS = AllDatRedReps.largest_ts>0;
            mnTS(j) = (sum(AllDatRedReps.largest_ts(AllDatRedReps.largest_ts>0))+...
                sum(AllDatRedReps.second_largest_ts(AllDatRedReps.second_largest_ts>0)))/sum(isTS);
        end
        dataResults.fracTSDatstd(iT) = std(frac);
        dataResults.meanNascentDatstd(iT) = std(mn);
        % dataResults.meanNascentTSDatstd(iT) = std(mnTS);
    end
end

% Adjust to account for partial transcripts.
nRNAmax = floor(ModelTS.fspOptions.bounds(nSp+3)+1);
PDO_partial_transcripts = zeros(nRNAmax,nRNAmax);
pBino = 0.5;
for i = 1:nRNAmax
    % PDO_partial_transcripts(1:i,i) = poisspdf(0:i-1,(i-1)/2);
    PDO_partial_transcripts(1:i,i) = binopdf(0:i-1,i-1,pBino);
end

% Calculate the likelihood function from the TS data.
TSLikelihood=0;
for iT=1:NT
    
    % Remove outliers from data and set to threshold.
    countsRNA = dataResults.TS_counts{iT};
    % countsRNA(countsRNA>outlierThresh) = outlierThresh;
    countsRNA = countsRNA(countsRNA<=outlierThresh);
   
    % Apply PDO to account for partial transcripts.
    P = PDO_partial_transcripts*modelResults.Psaved{iT};

    % Lump TS smaller than TSthresh as zero.
    P(1:TSthresh) = [sum(P(1:TSthresh));zeros(TSthresh-1,1)];

    % Bin outlier predictions for model
    % P = modelResults.Psaved{iT};

    P(outlierThresh+1) = sum(P(outlierThresh+1:end));
    P = P(1:outlierThresh+1);

    % Convolve two aleles to get sum of TS spots.
    Pconv = conv(P,P); 

    P = [Pconv(1:outlierThresh);sum(Pconv(outlierThresh+1:end))];

    modelResults.meanNascentMod(iT) = [TSthresh:length(P)-1]*P(TSthresh+1:end);
    modelResults.meanNascentTSMod(iT) = [TSthresh:length(P)-1]*P(TSthresh+1:end)/(1-sum(P(1:TSthresh)));
    x2 =  [TSthresh:length(P)-1].^2*P(TSthresh+1:end)/(1-sum(P(1:TSthresh)));
    modelResults.meanNascentTSModStd(iT) = sqrt(x2-modelResults.meanNascentTSMod(iT))./sqrt(dataResults.N(iT));
    modelResults.fracTSMod(iT) = sum(P(TSthresh+1:end));

    TSLikelihood = TSLikelihood + sum(log(P(countsRNA+1)));
    
end
end