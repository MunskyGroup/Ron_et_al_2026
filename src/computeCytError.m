function [TotError,ssaSoln] = computeCytError(x,indsPars,Model,saveBest,saveName, ...
    useGR,GRModel,parName,weights,makePlots,transformSSAsoln,...
    crnadat, ctime, creplica, crnadatGR, ctimeGR, creplicaGR)
arguments
    x
    indsPars
    Model
    saveBest = false
    saveName = ''
    useGR = false
    GRModel = []
    parName = 'bestx'
    weights = [1,1,1,1];
    makePlots = false
    transformSSAsoln = []
    crnadat = [34,4]
    ctime = 17
    creplica = 19
    crnadatGR = [23,24]
    ctimeGR = 5
    creplicaGR = 7
end
Model.parameters(indsPars,2) = num2cell(x);
ssaSoln = Model.solve(returnType='soln');

% Transformation to results of SSA.  For example, could combine two species
% into one.
if ~isempty(transformSSAsoln)
    ssaSoln.trajs = tensorprod(transformSSAsoln, ssaSoln.trajs, 2, 1);
end

[CytError,KSCyt] = compareDistPlots(ssaSoln,Model,6,2,[0:1:300],crnadat, ctime, creplica);
[NucError,KSNuc] = compareDistPlots(ssaSoln,Model,5,1,[0:1:300],crnadat, ctime, creplica);
TotError = CytError*weights(1) + NucError*weights(2);
if useGR
    [GRErrorCyt,KSGRCyt] = compareDistPlots(ssaSoln,GRModel,3,1,crnadatGR, ctimeGR, creplicaGR);
    [GRErrorNuc,KSGRNuc] = compareDistPlots(ssaSoln,GRModel,4,2,crnadatGR, ctimeGR, creplicaGR);
    TotError = CytError*weights(1) + NucError*weights(2) + GRErrorCyt*weights(3) + GRErrorNuc*weights(4);
end

if saveBest
    load(saveName,['minerr_',parName]);
    eval(['minerr=minerr_',parName,';']);
    if TotError<=minerr
        minerr = TotError
        eval(['minerr_',parName,'=minerr;']);
        eval([parName,' = x']);
        save(saveName,['minerr_',parName],parName,'-append')
        
        if makePlots
            close all
            [~,figs1] = makeCytDistPlots(ssaSoln,Model,601,5,1,[0:5:300],true,'cdf',crnadat, ctime, creplica,false);
            figure(figs1{2})
            [~,figs] = makeCytDistPlots(ssaSoln,Model,601,6,2,[0:5:300],true,'cdf',crnadat, ctime, creplica,false);
            figure(figs{2})
            % [~,figs3] = makeCytDistPlots(ssaSoln,GRModel,601,3,1,[0:40],true,'cdf',crnadatGR, ctimeGR, creplicaGR);
            % figure(figs3{2})
            % [~,figs4] = makeCytDistPlots(ssaSoln,GRModel,601,4,2,[0:40],true,'cdf',crnadatGR, ctimeGR, creplicaGR);
            % figure(figs4{2})
            drawnow
        end
    end
end

