function [distance,dists] = computeTotalObjFun(x,indsAllPars,SSAModel,TSLikelihood,...
    log10PriorMean,log10PriorStd,savePars,savedParsFile,parName,includeTS,...
    makePlots,defaultTSValue,weights_ABC,weights_CytNuc,transformSSAsoln,...
    crnadat, ctime, creplica, crnadatGR, ctimeGR, creplicaGR)
arguments
    x
    indsAllPars
    SSAModel
    TSLikelihood
    log10PriorMean =[]
    log10PriorStd =[]
    savePars = false
    savedParsFile = []
    parName = [];
    includeTS = true
    makePlots = false
    defaultTSValue = NaN;
    weights_ABC = [5000,500,10];
    weights_CytNuc = [1,1];
    transformSSAsoln = [];
    crnadat = [34,4]
    ctime = [17]
    creplica = [19]
    crnadatGR = [23,24]
    ctimeGR = 5
    creplicaGR = 7
end


[distance1,ssaSoln] = computeCytError(x,indsAllPars,SSAModel,false,'',...
    false,[],'',[weights_CytNuc,0,0],false,transformSSAsoln,...
    crnadat, ctime, creplica, crnadatGR, ctimeGR, creplicaGR);

if includeTS
    [distance2,TSmodelResults,TSdataResults] = TSLikelihood(x);
else
    distance2 = defaultTSValue*weights_ABC(2);
end
if ~isempty(log10PriorMean)&&~isempty(log10PriorStd)
    priorDistance = sum((log10(reshape(x,[numel(x),1]))-log10PriorMean(indsAllPars)).^2./(2*log10PriorStd(indsAllPars)));
else
    priorDistance=0;
end
distance = distance1/weights_ABC(1) - distance2/weights_ABC(2) + priorDistance/weights_ABC(3);

dists = [ -distance1/weights_ABC(1) , distance2/weights_ABC(2) , -priorDistance/weights_ABC(3)];

if savePars
    load(savedParsFile,['minerr_',parName]);
    eval(['minerr=minerr_',parName,';']);
    if distance<minerr
        minerr = distance
        dists
        eval(['minerr_',parName,'=minerr;']);
        eval([parName,'=x;'])
        save(savedParsFile, ['minerr_',parName], parName, '-append');

        if makePlots
            close all
            % Make Nuclear DUSP1 figure
            [~,figs1] = makeCytDistPlots(ssaSoln,SSAModel,601,5,1,[0:5:300],false,'cdf',crnadat, ctime, creplica,false);
            figure(figs1{2})
            % Make Cytoplasmic DUSP1 figure
            [~,figs] = makeCytDistPlots(ssaSoln,SSAModel,601,6,2,[0:5:300],false,'cdf',crnadat, ctime, creplica,false);
            figure(figs{2})
            % Plot TS

            % figure
            % X = [SSAModel_100.tSpan(2:end),SSAModel_100.tSpan(end:-1:2)];
            % Y = [TSmodelResults.meanNascentTSMod+TSmodelResults.meanNascentTSModStd,...
            %     TSmodelResults.meanNascentTSMod(end:-1:1)-TSmodelResults.meanNascentTSModStd(end:-1:1)];
            % fill(X,Y,[0.8 0.8 1]); hold on
            % plot(SSAModel_100.tSpan(2:end),TSmodelResults.meanNascentTSMod,'b','linewidth',3);
            % errorbar(SSAModel_100.tSpan(2:end),TSdataResults.meanNascentTSDat,TSdataResults.meanNascentTSDatstd)

            if includeTS
                dummy.tSpan = SSAModel.tSpan(2:end);
                figs = makeTSPlots(TSmodelResults,TSdataResults,dummy);
                figure(figs{2});
            end

            drawnow
        end
    end
end

end
