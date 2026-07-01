metricsFile = 'savedResults/Metrics.mat';
load(metricsFile)

%% Make MLE Comparison Table for GR
RedTable = GRModelComparison(:,[2,5,8,10]);
RedTable.DeltaLogPosterior = RedTable.DeltaLogPosterior;
RedTable.ApproxBF = -RedTable.DeltaBIC/2;
RedTable.ApproxBF(RedTable.ApproxBF==0) = 0;
RedTable.ApproxBF_Lapl = (GRModelComparison.LaplaceEvidence - max(GRModelComparison.LaplaceEvidence));

fid = fopen('ParameterTables/GR_ModelComparison.txt','w');
fprintf(fid,'\n%%Comparison of GR model fit errors:\n');
fprintf(fid,'\\begin{tabular}{|l||c|c|c|c|}');
fprintf(fid,'\\hline\n');
fprintf(fid,'& Import,$\\gamma_{\\rm nuc},\\gamma_{\\rm cyt}$ & Import,$\\gamma_{\\rm cyt}$ & Import,$\\gamma_{\\rm nuc}$ & Export,$\\gamma_{\\rm nuc},\\gamma_{\\rm cyt}$ \\\\ \\hline\\hline\n');

% rownames = {'Num.\ Parameters',...
%     '$\log_{10}(L/L_{\rm best})$',...
%     'AIC - AIC$_{\rm best}$', 'BIC - BIC$_{\rm best}$', '$\log_{10}(aBF_{\rm BIC})$', '$\log_{10}(aBF_{\rm Laplace})$'};
rownames = {'Num.\ Parameters',...
    '$\log(MAP/MAP_{\rm best})$',...
    'AIC - AIC$_{\rm best}$', 'BIC - BIC$_{\rm best}$', '$\log(aBF_{\rm BIC})$', '$\log(aBF_{\rm Laplace})$'};
for p = 1:5
    if p==1
        fprintf(fid,'%15s & %d & %d & %d & %d \\\\ \\hline\n', ...
            rownames{p}, RedTable{1,p}, RedTable{2,p}, RedTable{3,p}, RedTable{4,p});
    else
        % fprintf(fid,'%15s & %.2e & %.2e & %.2e & %.2e \\\\ \\hline\n', ...
        fprintf(fid,'%15s & %.0f & %.0f & %.0f & %.0f \\\\ \\hline\n', ...
            rownames{p}, RedTable{1,p}, RedTable{2,p}, RedTable{3,p}, RedTable{4,p});
    end
end
fprintf(fid,'\\end{tabular}');

%% Make MLE Comparison Table for DUSP1
RedTable = DUSP1NucModelComparison([1,2,3,5,4],[2,4,5,7,10,12]);
RedTable.LogPosterior_Nuc = (RedTable.LogPosterior_Nuc - max(RedTable.LogPosterior_Nuc));
RedTable.LogLikelihood_TS = (RedTable.LogLikelihood_TS - max(RedTable.LogLikelihood_TS));
RedTable.DeltaTotalLogPosterior = (RedTable.DeltaTotalLogPosterior - max(RedTable.DeltaTotalLogPosterior));

RedTable.DeltaBIC = RedTable.DeltaBIC - min(RedTable.DeltaBIC);
RedTable.DeltaAIC = RedTable.DeltaAIC - min(RedTable.DeltaAIC);

RedTable.ApproxBF = -RedTable.DeltaBIC/2;
RedTable.ApproxBF(RedTable.ApproxBF==0) = 0;

fid = fopen('ParameterTables/NucDUSP1_ModelComparison.txt','w');
fprintf(fid,'\n%%Comparison of models for DUSP1 TS and Nuclear mRNA models:\n');
fprintf(fid,'\\begin{tabular}{|l||c|c|c|c|c|c|}');
fprintf(fid,'\\hline\n');
% fprintf(fid,'& Import,$\\gamma_{\\rm nuc},\\gamma_{\\rm cyt}$ & Import,$\\gamma_{\\rm cyt}$ & Import,$\\gamma_{\\rm nuc}$ & Export,$\\gamma_{\\rm nuc},\\gamma_{\\rm cyt}$ \\\\ \\hline\\hline\n');
fprintf(fid,'& $K_{OFF}$ & $K_{ON}$ & $K_{R}$& $K_{OFF},K_{ON}$ & $K_{OFF},K_{ON},K_{R}$  \\\\ \\hline\\hline\n');

rownames = {'Num.\ Parameters',...
    '$\log(MAP/MAP_{\rm best})_{Nuc}$',...
    '$\log(MAP/MAP_{\rm best})_{TS}$',...
    '$\log(MAP/MAP_{\rm best})_{Tot}$',...
    'AIC - AIC$_{\rm best}$', 'BIC - BIC$_{\rm best}$', '$\log(BF_{\rm approx})$'};
for p = 1:7
    if p==1
        fprintf(fid,'%15s & %d & %d & %d & %d & %d \\\\ \\hline\n', ...
            rownames{p}, RedTable{1,p}, RedTable{2,p}, RedTable{3,p}, RedTable{4,p}, RedTable{5,p});
    else
        % fprintf(fid,'%15s & %.2e & %.2e & %.2e & %.2e & %.2e \\\\ \\hline\n', ...
        fprintf(fid,'%15s & %.0f & %.0f & %.0f & %.0f & %.0f \\\\ \\hline\n', ...
            rownames{p}, RedTable{1,p}, RedTable{2,p}, RedTable{3,p}, RedTable{4,p}, RedTable{5,p});
    end
end
fprintf(fid,'\\end{tabular}');


%% Make Fit and Prediciton Error Tables for Full Models
ModelNamesFiles = {'Pars_first_order_decay','FinalModelNoTimeVariation','FinalModelNoSaturation','FinalModel_allButGR','Pars_GR_DUSP1_TTP_Model_All'};
ModelDisplayNames = {'Linear \& Constant','Nonlinear \& Constant','Linear \& Time-Varying','Nonlinear \& Time-Varying','Detailed Mechanistic'};

T = struct();
T(1).DataName = 'Nuc/Cyt, 100 nM Dex (Fit)';
T(2).DataName = 'TS, 100 nM Dex (Fit)';
T(3).DataName = 'Nuc/Cyt, 0 min Tpl (Fit)';

% Dex at 0.3, 1, 10 nM
T(4).DataName = 'Nuc/Cyt, 0.3 nM Dex (Predict)';
T(5).DataName = 'Nuc/Cyt, 1.0 nM Dex (Predict)';
T(6).DataName = 'Nuc/Cyt, 10 nM Dex (Predict)';

T(7).DataName = 'TS, 0.3 nM Dex (Predict)';
T(8).DataName = 'TS, 1.0 nM Dex (Predict)';
T(9).DataName = 'TS, 10 nM Dex (Predict)';

% Tpl at 0, 20, 50, 75, 150, 180
T(10).DataName = 'Nuc/Cyt, 20 min Tpl (Predict)';
T(11).DataName = 'Nuc/Cyt, 75 min Tpl (Predict)';
T(12).DataName = 'Nuc/Cyt, 150 min Tpl (Predict)';
T(13).DataName = 'Nuc/Cyt, 180 min Tpl (Predict)';
T(14).DataName = 'Total Fit Error';
T(15).DataName = 'Total Prediction Error';

for iM = 1:5
    % Fitting Data
    T(1).(ModelNamesFiles{iM}) = eval(['JNucCyt_Predict_',ModelNamesFiles{iM},'_Dex_100']); %'Nuc/Cyt, 100 nM Dex (Fit)';
    T(2).(ModelNamesFiles{iM}) = eval(['JTS_Predict_',ModelNamesFiles{iM},'_Dex_100']); %'TS, 100 nM Dex (Fit)';
    T(3).(ModelNamesFiles{iM}) = eval(['JNucCyt_Predict_',ModelNamesFiles{iM},'_Tpl_0min']); %'Nuc/Cyt, 0 min Tpl (Fit)';

    % Dex at 0.3, 1, 10 nM
    T(4).(ModelNamesFiles{iM}) = eval(['JNucCyt_Predict_',ModelNamesFiles{iM},'_Dex_0p3']); %'Nuc/Cyt, 0.3 nM Dex (Predict)';
    T(5).(ModelNamesFiles{iM}) = eval(['JNucCyt_Predict_',ModelNamesFiles{iM},'_Dex_1']); %'Nuc/Cyt, 1.0 nM Dex (Predict)';
    T(6).(ModelNamesFiles{iM}) = eval(['JNucCyt_Predict_',ModelNamesFiles{iM},'_Dex_10']); %'Nuc/Cyt, 10 nM Dex (Predict)';

    T(7).(ModelNamesFiles{iM}) = eval(['JTS_Predict_',ModelNamesFiles{iM},'_Dex_0p3']); %'TS, 0.3 nM Dex (Predict)';
    T(8).(ModelNamesFiles{iM}) = eval(['JTS_Predict_',ModelNamesFiles{iM},'_Dex_1']); %'TS, 1.0 nM Dex (Predict)';
    T(9).(ModelNamesFiles{iM}) = eval(['JTS_Predict_',ModelNamesFiles{iM},'_Dex_10']); %'TS, 10 nM Dex (Predict)';

    % Tpl at 0, 20, 50, 75, 150, 180
    T(10).(ModelNamesFiles{iM}) = eval(['JNucCyt_Predict_',ModelNamesFiles{iM},'_Tpl_20min']); %'Nuc/Cyt, 20 min Tpl (Predict)';
    T(11).(ModelNamesFiles{iM}) = eval(['JNucCyt_Predict_',ModelNamesFiles{iM},'_Tpl_75min']); %'Nuc/Cyt, 75 min Tpl (Predict)';
    T(12).(ModelNamesFiles{iM}) = eval(['JNucCyt_Predict_',ModelNamesFiles{iM},'_Tpl_150min']); %'Nuc/Cyt, 150 min Tpl (Predict)';
    T(13).(ModelNamesFiles{iM}) = eval(['JNucCyt_Predict_',ModelNamesFiles{iM},'_Tpl_180min']); %'Nuc/Cyt, 180 min Tpl (Predict)';
end

for p = 1:13
    best = max([T(p).Pars_first_order_decay, T(p).FinalModelNoTimeVariation, ...
       T(p).FinalModelNoSaturation, T(p).FinalModel_allButGR, T(p).Pars_GR_DUSP1_TTP_Model_All]);
    T(p).Pars_first_order_decay = T(p).Pars_first_order_decay - best;
    T(p).FinalModelNoTimeVariation = T(p).FinalModelNoTimeVariation - best;
    T(p).FinalModelNoSaturation = T(p).FinalModelNoSaturation - best;
    T(p).FinalModel_allButGR = T(p).FinalModel_allButGR - best;
    T(p).Pars_GR_DUSP1_TTP_Model_All = T(p).Pars_GR_DUSP1_TTP_Model_All - best;
end

for iM = 1:5
    T(14).(ModelNamesFiles{iM}) = sum([T(1:3).(ModelNamesFiles{iM})]);
    T(15).(ModelNamesFiles{iM}) = sum([T(4:13).(ModelNamesFiles{iM})]);
end

% Write to LaTeX
fid = fopen('ParameterTables/FitAndPredictionErrors.txt','w');
fprintf(fid,'\n%%Comparison of model fit and prediction errors:\n');
fprintf(fid,'%%Data Set & Linear & Constant & Nonlinear & Constant & Linear & Time-Varying & Nonlinear & Time-Varying & Detailed Mechanistic \\hline \\\\ \n');
fprintf(fid,'\\begin{tabular}{|l||l|l|l|l|l|}\n')
fprintf(fid,'\\hline\n')
fprintf(fid,'\\textbf{Data Set} & \\makecell{\\textbf{Linear}\\\\ \\textbf{Constant}} & \\makecell{\\textbf{Nonlinear}\\\\ \\textbf{Constant}} & \\makecell{\\textbf{Linear}\\\\ \\textbf{Time-Varying}} & \\makecell{\\textbf{Nonlinear}\\\\ \\textbf{Time-Varying}} & \\makecell{\\textbf{Detailed}\\\\ \\textbf{Mechanistic}} \\\\ \\hline \\hline\n')
for p = 1:3 
   fprintf(fid,'%15s & %.0f & %.0f & %.0f & %.0f & %.0f \\\\ \\hline\n', ...
       T(p).DataName, T(p).Pars_first_order_decay, T(p).FinalModelNoTimeVariation, ...
       T(p).FinalModelNoSaturation, T(p).FinalModel_allButGR, T(p).Pars_GR_DUSP1_TTP_Model_All);
end
p=14;
fprintf(fid,'{\\bf %15s} & %.0f & %.0f & %.0f & %.0f & %.0f \\\\ \\hline\\hline\n', ...
    T(p).DataName, T(p).Pars_first_order_decay, T(p).FinalModelNoTimeVariation, ...
    T(p).FinalModelNoSaturation, T(p).FinalModel_allButGR, T(p).Pars_GR_DUSP1_TTP_Model_All);

for p = 4:13
    fprintf(fid,'%15s & %.0f & %.0f & %.0f & %.0f & %.0f \\\\ \\hline\n', ...
        T(p).DataName, T(p).Pars_first_order_decay, T(p).FinalModelNoTimeVariation, ...
        T(p).FinalModelNoSaturation, T(p).FinalModel_allButGR, T(p).Pars_GR_DUSP1_TTP_Model_All);
end
p=15;
fprintf(fid,'{\\bf %15s} & %.0f & %.0f & %.0f & %.0f & %.0f \\\\ \\hline\n', ...
    T(p).DataName, T(p).Pars_first_order_decay, T(p).FinalModelNoTimeVariation, ...
    T(p).FinalModelNoSaturation, T(p).FinalModel_allButGR, T(p).Pars_GR_DUSP1_TTP_Model_All);

fprintf(fid,'\\end{tabular}')
fclose(fid);

