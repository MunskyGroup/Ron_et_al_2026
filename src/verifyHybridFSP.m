%% VerifyHybridFSP
load savedModels/KON.mat
KON = KON.formPropensitiesGeneral('KON');

%% Add negative initial time for SS approx.
KON.initialTime = -2000;
KON.tSpan = [KON.tSpan];
KON.solutionScheme = 'fsp';
[~,~,KON] = KON.solve;

%%
try
    KON.Solutions = rmfield(KON.Solutions,'trajs');
catch
end
KON.ssaOptions.Nsims = 1000000;
KON.solutionScheme = 'ssa';
[~,~,KON] = KON.solve;

%% Call function to make verification plots.
KON = KON.verifyFSPandSSA(speciesNames={'rna'});

set(gca,'XLim',[0,250]);
f = gcf;
legend('1,000,000 Exact SSA Samples','Hybrid ODE-FSP');
f.Children(3).Title.String = '';
saveas(f,'savedFigures/FinalFigures/SupplementalFigures/VerifyFSP.fig');
f.Children(3).YScale = 'log';
saveas(f,'savedFigures/FinalFigures/SupplementalFigures/VerifyFSPLog.fig');
%%
for i = 1:100
    KON.Solutions = rmfield(KON.Solutions,'trajs');
    KON.ssaOptions.Nsims = 1000;
    KON.solutionScheme = 'ssa';
    [~,~,KON] = KON.solve;

    P = double(KON.Solutions.fsp{end}.p.sumOver([1,2]).data);
    D = squeeze(KON.Solutions.trajs(5,end,:));
    [~,p(i)] = kstest(D,"CDF",[(0:length(P)-1)',cumsum(P)]);
end
median(p)