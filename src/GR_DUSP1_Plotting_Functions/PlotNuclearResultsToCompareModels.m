close all

dataFileDusp1 = 'Data/DUSP1_SSITcellresults_Final_Sep18.csv';
modelLibrary = 'savedParameters/GRDusp1ModelLibrary';
savedParsFile ='savedParameters/Parameters_Final';

% Load data for Transcription Sites.
allData = readtable(dataFileDusp1);
allData = allData(allData.cyto_area_px>=12593&allData.cyto_area_px<=17685&isnan(allData.time_tpl),:);

dex = 100;
errorValues = [];

load(savedParsFile,'GRpars');

savedNucParSet = {'Pars_KOFF','Pars_KON','Pars_KR','Pars_KON_KOFF_KR',...
    'Pars_KON_KOFF','Pars_KON_KOFF_KR_zeroOFF','Pars_KON_KOFF_zeroOFF','Pars_KON_zeroOFF'};

for iM = 1:3
    [Mod,log10PriorMean,log10PriorStd] = ...
        dusp1ModelLibrary_Final('ModelDUSP1_100nM',false);

    Mod.parameters(5:12,2) = num2cell(GRpars);

    load(savedParsFile,savedNucParSet{iM})
    eval(['DUSP1parsAlt=',savedNucParSet{iM},';']);

    switch iM
        case 1 % KOFF
            Mod.fittingOptions.modelVarsToFit = [1,2,3,4,19,23,25];
            zeroPars = [24,26];
        case 2 % KON
            Mod.fittingOptions.modelVarsToFit = [1,2,3,4,19,23,24];
            zeroPars = [25,26];
        case 3 % KR
            Mod.fittingOptions.modelVarsToFit = [1,2,3,4,19,23,26];
            zeroPars = [24,25];
    end

    % case 3
    %     Mod = load('savedModels/KR.mat'); Mod = Mod.KR;
    % case 2
    %     Mod = load('savedModels/KON.mat'); Mod = Mod.KON;
    % case 1
    %       Mod = load('savedModels/KOFF.mat'); Mod = Mod.KOFF;
    % % KOFF seems to be missing, so I will just load the parameters
    %   % directly.
    %   Mod = load('savedModels/KON.mat'); Mod = Mod.KON;
    %   Mod.fittingOptions.modelVarsToFit = [1,2,3,4,19,23,25];
    %   zeroPars = [24,26];
    %   Mod.parameters(zeroPars,2) = num2cell(0*zeroPars);
    %   savedParsFile ='savedParameters/Parameters_Final';
    %   load(savedParsFile,'Pars_KOFF')
    %   Mod.parameters(Mod.fittingOptions.modelVarsToFit,2) = num2cell(Pars_KOFF);
    % case 4
    %     Mod = load('savedModels/KON_KOFF_KR.mat'); Mod = Mod.KON_KOFF_KR;

    % end

    Mod.parameters(13,2) = {dex};
    % Update parameters from previous saved file.
    Mod.parameters(Mod.fittingOptions.modelVarsToFit,2) = num2cell(DUSP1parsAlt);

    % Set removed parameters to zero.
    Mod.parameters(zeroPars,2) = num2cell(zeros(length(zeroPars),1));

    Mod = Mod.solve;

    Mod.fittingOptions.logPrior = [];

    [errorValues(iM,1),~,Results{iM}] = Mod.computeLikelihood;
    for iT = 1:length(Mod.Solutions.fsp)
        P = squeeze(sum(double(Mod.Solutions.fsp{iT}.p.data),[1,2]));
        mu(iT) = [0:length(P)-1]*P;
        mu2(iT) = [0:length(P)-1].^2*P;
        var(iT) = mu2(iT) - mu(iT)^2;
    end

    % Update model for TS analysis
    ModelTS = Mod;
    ModelTS.tSpan = unique(allData.time)';
    ModelTS.fittingOptions.logPrior = []; % Remove parameters prior since this is already included in the nuclear model.

    % tsObj = @(x)computeTSlikelihoodSSA(x,Mod,allData,50,true);
    tsObj = @(x)computeTSlikelihood(x,ModelTS,allData,50,true);

    [errorValues(iM,2),modelResults,dataResults] = tsObj([Mod.parameters{Mod.fittingOptions.modelVarsToFit,2}]);

    J = dataResults.N~=0;

    figure(1)
    subplot(2,2,1)
    plot(Mod.tSpan,modelResults.fracTSMod,'LineWidth',3); hold on;
    set(gca,'yLim',[0,1])

    subplot(2,2,2)
    plot(Mod.tSpan,modelResults.meanNascentTSMod,'LineWidth',3); hold on;
    set(gca,'yLim',[0,15])

    subplot(2,2,3)
    plot(Mod.tSpan(J),modelResults.meanNascentTSModStd(J).*sqrt(dataResults.N(J)),'LineWidth',3); hold on;

    subplot(2,2,4)
    plot(Mod.tSpan(J),(modelResults.meanNascentTSModStd(J).*sqrt(dataResults.N(J))).^2./modelResults.meanNascentTSMod(J),'LineWidth',3); hold on;

    figure(2)
    subplot(1,2,1)
    plot(Mod.tSpan,mu,'LineWidth',3); hold on;
    % set(gca,'yLim',[0,1])

    subplot(1,2,2)
    plot(Mod.tSpan,sqrt(var),'LineWidth',3); hold on;
    % set(gca,'yLim',[0,15])


end

figure(1)
subplot(2,2,1)
title('Fraction of Cells with TS')
xlabel('Time (min)')
ylabel('Fraction Active')
plot(Mod.tSpan,dataResults.fracTSDat,'s'); hold on;
set(gca,'yLim',[0,1])
legend({'KOFF','KON','KR','ALL'},"Location","northeast")

subplot(2,2,2)
title('Mean Nascent RNA per Active TS')
xlabel('Time (min)')
ylabel('Mean RNA per ATS')
plot(Mod.tSpan,dataResults.meanNascentTSDat,'s'); hold on;
set(gca,'yLim',[0,15])

subplot(2,2,3)
title('STD Nascent RNA per Active TS')
xlabel('Time (min)')
ylabel('STD RNA per ATS')
plot(Mod.tSpan,dataResults.meanNascentDatstd.*sqrt(dataResults.N),'s'); hold on;
set(gca,'yLim',[0,15])

subplot(2,2,4)
title('Fano Factor Nasc. per Active TS')
xlabel('Time (min)')
ylabel('Fano Factor')
plot(Mod.tSpan,(dataResults.meanNascentDatstd.*sqrt(dataResults.N)).^2./dataResults.meanNascentTSDat,'s'); hold on;
set(gca,'yLim',[0,30])

figure(2)
subplot(1,2,1)
title('Mean Nuclear RNA')
plot(Mod.tSpan,Mod.dataSet.mean,'s')
set(gca,'yLim',[0,60])

subplot(1,2,2)
title('STD Nuclear RNA')
plot(Mod.tSpan,sqrt(Mod.dataSet.var),'s')
set(gca,'yLim',[0,35])


%%
figure(3)
plot(Results{1}.DataLoadingAndFittingTabOutputs.V_LogLk-Results{1}.DataLoadingAndFittingTabOutputs.V_LogLk); hold on
plot(Results{2}.DataLoadingAndFittingTabOutputs.V_LogLk-Results{1}.DataLoadingAndFittingTabOutputs.V_LogLk)
plot(Results{3}.DataLoadingAndFittingTabOutputs.V_LogLk-Results{1}.DataLoadingAndFittingTabOutputs.V_LogLk)
% plot(Results{4}.DataLoadingAndFittingTabOutputs.V_LogLk-Results{1}.DataLoadingAndFittingTabOutputs.V_LogLk)

figure(4)
plot(Results{1}.DataLoadingAndFittingTabOutputs.numCells)