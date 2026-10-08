%% Predictions Made by Full Mechanistic Model
addpath(genpath('src'));

%% Load Mechanistic Model and Update to Latest Parameters.
MechModel = dusp1ModelLibrary_Final('GR_DUSP1_TTP_Model_2',false);

savedParsFile ='savedParameters/Parameters_Final';
indsFreePars =  [1,2,3,4,14,15,19,20,21,22,24,25,26,27:36];
MechModel.fittingOptions.modelVarsToFit = indsFreePars;
load(savedParsFile,'Pars_GR_DUSP1_TTP_Model_All_2','GRpars');
valsFreePars = Pars_GR_DUSP1_TTP_Model_All_2;
zeroPars = [20,21,22,23];
includeTS = true;
transformSSAsoln = zeros(6,13);
transformSSAsoln(1:6,1:6) = eye(6);
transformSSAsoln(6,12) = 1;
MechModel.parameters(5:12,2) = num2cell(GRpars);
MechModel.parameters(indsFreePars,2) = num2cell(valsFreePars);

%% Modify to remove protection of protein degradation.
% MechModel = MechModel.addReaction(struct(...
%         'propensity','degProt*TTP_rCyt_complex'  ,...
%         'stoichiometry',{{'TTP_rCyt_complex',-1}}));
% 
% MechModel = MechModel.addReaction(struct(...
%         'propensity','degProt*rCyt_complex'  ,...
%         'stoichiometry',{{'rCyt_complex',-1}}));

MechModel = MechModel.changeParameter({'degProt',0.005});

MechModel = MechModel.changeParameter({'TTP_eta',1.0});

MechModel = MechModel.addParameter({'tPuro',200});
MechModel.propensityFunctions{5} = 'kg1*(t<tPuro)';
MechModel.propensityFunctions{15} = 'ktranslate*TTP_rCyt*(t<tPuro)';

%% Solve SSA 
figure(1); clf;
tParray = [200,0,75];
cols = {'k','b','r'};

load('savedParameters/MHResults_Pars_GR_DUSP1_TTP_Model_All_2')

for j = 1:1
    iSamp = 5000+randi([0,5000]);
    MechModel.parameters(indsFreePars,2) = num2cell(exp(MHResultsAll.mhSamples(iSamp,:)));
    
    % Change Deg parameter
    MechModel = MechModel.changeParameter({'degProt',0.005});

    for i = 1:3

        MechModel.ssaOptions.Nsims = 1e4;
        MechModel = MechModel.changeParameter({'tPuro',tParray(i)});
        MechModel = MechModel.solve(solver='ssa');

        % Make Plots of Various Species vs time
        t = MechModel.tSpan(2:end);
        s = MechModel.Solutions.trajs(:,2:end,:);

        nG = squeeze(s(4,:,:));
        cG = squeeze(s(3,:,:));
        nD = squeeze(s(5,:,:));
        cD = squeeze(s(6,:,:)+s(12,:,:));
        cDu = squeeze(s(6,:,:));
        cDb = squeeze(s(12,:,:));
        % cD = squeeze(s(6,:,:));
        nT=squeeze(s(9,:,:));
        cT=squeeze(s(10,:,:)+s(13,:,:));
        cTu=squeeze(s(10,:,:));
        cTb=squeeze(s(13,:,:));
        % cT=squeeze(s(10,:,:));
        pT=squeeze(s(11,:,:)+s(12,:,:)+s(13,:,:));
        % pT=squeeze(s(11,:,:)+s(13,:,:));
        % pT=squeeze(s(11,:,:));

        figure(1);

        subplot(4,3,1)
        plot(t,mean(nG,2),cols{i}); hold on
        ylabel('Nuc GR')

        subplot(4,3,2)
        plot(t,mean(cG,2),cols{i}); hold on
        ylabel('Cyt GR')

        subplot(4,3,3)
        plot(t,mean(cG+nG,2),cols{i}); hold on
        ylabel('Total GR')

        subplot(4,3,4)
        plot(t,mean(nD,2),cols{i}); hold on
        ylabel('Nuc DUSP1 RNA')

        subplot(4,3,5)
        plot(t,mean(cD,2),cols{i}); hold on
        ylabel('Cyt DUSP1 RNA')

        subplot(4,3,6)
        plot(t,mean(cDu,2)./mean(cD,2),cols{i}); hold on
        ylabel('Frac Unbound')

        subplot(4,3,7)
        plot(t,mean(nT,2),cols{i}); hold on
        ylabel('Nuc TTP RNA')

        subplot(4,3,8)
        plot(t,mean(cT,2),cols{i}); hold on
        ylabel('Cyt TTP RNA')

        subplot(4,3,9)
        plot(t,mean(cTu,2)./mean(cT,2),cols{i}); hold on
        ylabel('Frac Unbound')

        subplot(4,3,10)
        plot(t,mean(pT,2),cols{i}); hold on
        ylabel('TTP Protein')

    end
end
%% Correlation plots.
figure(2); clf;
sps = {cG,nG,nD,cD,nT,cT,pT};
legs = {'GR_c','GR_n','D_n','D_c','T_n','T_c','T_p'};
iT = 12;
for i =1:6
    for j = i+1:7
        subplot(6,6,(i-1)*6+j-1)
        rplot(sps{j}(iT,:),sps{i}(iT,:),'o')
        if j==i+1
            ylabel(legs{i})
            xlabel(legs{j})
        end
    end
end

%% Covarance
MAT = [cG(iT,:);nG(iT,:);nD(iT,:);cD(iT,:);nT(iT,:);cT(iT,:);pT(iT,:)];
COV = cov(MAT');
for i =1:6
    for j = i+1:7
        disp(['COV(',legs{i},',',legs{j},')=',num2str(COV(i,j)/sqrt(COV(i,i)*COV(j,j)))])
    end
end

%%
function rplot(x,y,m)
plot(x-0.5+rand(size(x)),y-0.5+rand(size(x)),m)
end