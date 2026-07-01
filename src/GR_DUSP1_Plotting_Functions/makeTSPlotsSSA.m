function figHandles = makeTSPlotsSSA(modelResults,dataResults,ModelTS,TSthresh,dPlotType)
arguments
    modelResults
    dataResults
    ModelTS
    TSthresh = 4;
    dPlotType = 'pdf'
end
inds = find(~isnan(dataResults.fracTSDat));

figHandles{1} = figure('Visible','on');
NT = length(modelResults.meanNascentMod);
NTd = length(inds);
for iT = 1:NTd
    P=modelResults.Psaved{inds(iT)};
    subplot(ceil(NTd/3),3,iT)
    switch dPlotType 
        case 'cdf'
            P(1:TSthresh) = [sum(P(1:TSthresh)),zeros(1,TSthresh-1)];  
            stairs([0:length(P)-1],cumsum(P),'linewidth',3)
            set(gca,'fontsize',16,'ylim',[0,1],'xlim',[0,50])
        case 'pdf'
            P(1:TSthresh) = sum(P(1:TSthresh))/4;
            stairs([0:length(P)-1],(P),'linewidth',3)
            set(gca,'fontsize',16,'ylim',[0,0.2],'xlim',[0,50])
    end
    if iT>9
        xlabel('Number nascent RNA')
    end
    if mod(iT-1,3)==0
        ylabel('Probability')
    end

end

for iT = 1:NTd
    subplot(ceil(NTd/3),3,iT); hold on
    histogram(dataResults.TS_counts{inds(iT)},[0,TSthresh:52],'Normalization',dPlotType,'DisplayStyle','bar','LineWidth',1);
    % histogram(dataResults.TS_counts{iT},[0,TSthresh:52],'Normalization','cdf','DisplayStyle','stairs','LineWidth',3);
    % nCount(iT) = length(dataResults.TS_counts{iT});
    % semCount(iT) = std(dataResults.TS_counts{iT})/sqrt(nCount(iT)-1);
end


figHandles{2} = figure('Visible','on');
subplot(2,1,1)
T = [ModelTS.tSpan,ModelTS.tSpan(end:-1:1)];
Y = [modelResults.meanNascentTSMod-modelResults.meanNascentTSModStd,...
    modelResults.meanNascentTSMod(end:-1:1)+modelResults.meanNascentTSModStd(end:-1:1)];
Y(isnan(Y)) = 0;
fill(T(isfinite(Y)),Y(isfinite(Y)),'b');hold on
plot(ModelTS.tSpan,modelResults.meanNascentTSMod,'linewidth',3);
% set(gca,'fontsize',16,'ylim',[0,5])
ylabel('Mean Nascent RNA')

subplot(2,1,2)
% plot(ModelTS.tSpan,modelResults.fracTSMod,'linewidth',3);
stdExp = 1./sqrt(dataResults.N).*sqrt(modelResults.fracTSMod'.*(1-modelResults.fracTSMod'));
Y = [modelResults.fracTSMod'-stdExp,...
    modelResults.fracTSMod(end:-1:1)'+stdExp(end:-1:1)];
fill(T(isfinite(Y)),Y(isfinite(Y)),'b');hold on
plot(ModelTS.tSpan,modelResults.fracTSMod,'linewidth',3);
set(gca,'fontsize',16,'ylim',[0,0.6])
ylabel('Fraction with TS')
xlabel('time (min)')
%
%    STEP 4.C.2. -- Add data to plots

subplot(2,1,1)
hold on
% plot(ModelTS.tSpan,dataResults.meanNascentDat,'s','markersize',16,'linewidth',3,'MarkerFaceColor','k');
% errorbar(ModelTS.tSpan,dataResults.meanNascentDat,dataResults.meanNascentDatstd,'s','markersize',16,'linewidth',3,'MarkerFaceColor','k');
errorbar(ModelTS.tSpan,dataResults.meanNascentTSDat,dataResults.meanNascentTSDatstd,'s','markersize',16,'linewidth',3,'MarkerFaceColor','k');

subplot(2,1,2)
hold on
% stdFrac = sqrt(nCount.*dataResults.fracTSDat.*(1-dataResults.fracTSDat))./nCount;
% plot(ModelTS.tSpan,dataResults.fracTSDat,'s','markersize',16,'linewidth',3,'MarkerFaceColor','k');
errorbar(ModelTS.tSpan,dataResults.fracTSDat,dataResults.fracTSDatstd,'s','markersize',16,'linewidth',3,'MarkerFaceColor','k');
end
