function Wrapper_For_Final_Model_Comparisons(opts)
arguments
    opts
end
%% Define models and labels
load("savedModels/FinalModel_allButGR_NoPrior")
load("savedModels/FinalModelNoSaturation")
load("savedModels/FinalModelNoTimeVariation")
load("savedModels/First_order_decay")
load("savedModels/GR_DUSP1_TTP_Model_All_2")

modelSpecs = struct([]);

modelSpecs(1).Model = FinalModel_allButGR_NoPrior;
modelSpecs(1).Label = "T.V. + Sat.";
modelSpecs(1).transformSSAsoln = [];

modelSpecs(2).Model = FinalModelNoSaturation;
modelSpecs(2).Label = "T.V. Only";
modelSpecs(2).transformSSAsoln = [];

modelSpecs(3).Model = FinalModelNoTimeVariation;
modelSpecs(3).Label = "Sat. Only";
modelSpecs(3).transformSSAsoln = [];

modelSpecs(4).Model = First_order_decay;
modelSpecs(4).Label = "First Order Decay";
modelSpecs(4).transformSSAsoln = [];

modelSpecs(5).Model = GR_DUSP1_TTP_Model_All_2;
modelSpecs(5).Label = "Mechanistic TTP";
modelSpecs(5).transformSSAsoln = zeros(6,13);
modelSpecs(5).transformSSAsoln(1:6,1:6) = eye(6);
modelSpecs(5).transformSSAsoln(6,12) = 1;

%% Species to analyze

for iComp = 1:2
    switch iComp
        case 1
            % Cytoplasmic DUSP1 mRNA
            rIdxModel = 6;
            speciesLabel = "mRNA_{cyt}";
            rIdxData = 2;
            saveFileName = "savedResults/DUSP1_cyt_final_models_fano.csv";
        case 2
            % Nuclear DUSP1 mRNA
            rIdxModel = 5;
            speciesLabel = "mRNA_{nuc}";
            rIdxData = 1;
            saveFileName = "savedResults/DUSP1_nuc_final_models_fano.csv";
    end


    %% Experimental data set
    % Use the corresponding cytoplasmic dataSet.
    % For example:
    expDataSet = FinalModel_allButGR_NoPrior.dataSet;

    %% Run analysis and make plots

    [T_cyt_fano, figs{iComp}] = fanosFromSSATrajs( ...
        modelSpecs, ...
        expDataSet, ...
        rIdxModel, ...
        speciesLabel, ...
        rIdxData, ...
        TimePoint = 75, ...
        ExcludeTimesBelow = 0, ...
        MakePlots = true, ...
        SaveFigures = false, ...
        OutDir = "savedFigures/Moments", ...
        FilePrefix = "DUSP1_nuc_final_models" , ...
        showErrorBars = false);

    disp(T_cyt_fano)

    writetable(T_cyt_fano, saveFileName);

end

%%
F1 = figs{2}.meanVsTime; figure(F1)
set(F1,"Position",[428   677   264   204]);
F1.Children(1).Visible = 'off';
F1.Children(2).FontSize = 10;
F1.Children(2).FontWeight = 'bold';
F1.Children(2).XLabel.Visible = 'off';
F1.Children(2).YLabel.Visible = 'off';
F1.Children(2).XLim = [0,185];
F1.Children(2).YLim = [15,60];
F1.Children(2).XTick = [0:30:180];
F1.Children(2).YTick = [15:15:60];
F1.Children(2).XTickLabelRotation = 0;

F2 = figs{2}.varVsTime; figure(F2)
set(F2,"Position",[428   677   264   204]);
F2.Children(1).Visible = 'off';
F2.Children(2).FontSize = 10;
F2.Children(2).FontWeight = 'bold';
F2.Children(2).XLabel.Visible = 'off';
F2.Children(2).YLabel.Visible = 'off';
F2.Children(2).XLim = [0,185];
F2.Children(2).YLim = [0,1600];
F2.Children(2).XTick = [0:30:180];
F2.Children(2).YTick = [0:400:1600];
F2.Children(2).XTickLabelRotation = 0;

F3 = figs{1}.meanVsTime; figure(F3)
set(F3,"Position",[428   677   264   204]);
F3.Children(1).Visible = 'off';
F3.Children(2).FontSize = 10;
F3.Children(2).FontWeight = 'bold';
F3.Children(2).XLabel.Visible = 'off';
F3.Children(2).YLabel.Visible = 'off';
F3.Children(2).XLim = [0,185];
F3.Children(2).YLim = [30,130];
F3.Children(2).XTick = [0:30:180];
F3.Children(2).YTick = [30:30:120];
F3.Children(2).XTickLabelRotation = 0;

F4 = figs{1}.varVsTime; figure(F4)
set(F4,"Position",[428   677   264   204]);
F4.Children(1).Visible = 'off';
F4.Children(2).FontSize = 10;
F4.Children(2).FontWeight = 'bold';
F4.Children(2).XLabel.Visible = 'off';
F4.Children(2).YLabel.Visible = 'off';
F4.Children(2).XLim = [0,185];
F4.Children(2).YLim = [0,4000];
F4.Children(2).XTick = [0:30:180];
F4.Children(2).YTick = [0:1000:4000];
F4.Children(2).XTickLabelRotation = 0;

for i = 2:6
    F1.Children(2).Children(i).Marker = "none";
    F2.Children(2).Children(i).Marker = "none";
    F3.Children(2).Children(i).Marker = "none";
    F4.Children(2).Children(i).Marker = "none";
end

saveas(F1,'savedFigures/FinalFigures/MainFigures/Fig3DTL_MeanNuc_ModelComparison.fig','fig');
saveas(F2,'savedFigures/FinalFigures/MainFigures/Fig3DTR_VarNuc_ModelComparison.fig','fig');
saveas(F3,'savedFigures/FinalFigures/MainFigures/Fig3DBL_MeanCyt_ModelComparison.fig','fig');
saveas(F4,'savedFigures/FinalFigures/MainFigures/Fig3DBR_VarCyt_ModelComparison.fig','fig');

%% Display the means and fano factors to the user.
Data = T_cyt_fano(strcmp(T_cyt_fano.Type,"Data"),:);

disp(['FANO Factor = ',num2str(Data.Fano(1)),' pm ',num2str(Data.FanoStd(1)),' at t = ',num2str(0),' min'])
disp(['Mean at same time = ',num2str(Data.Mean(1)),' pm ',num2str(Data.MeanStd(1))])

disp(['FANO Factor = ',num2str(Data.Fano(9)),' pm ',num2str(Data.FanoStd(9)),' at t = ',num2str(90),' min'])
disp(['Mean at same time = ',num2str(Data.Mean(9)),' pm ',num2str(Data.MeanStd(9))])