function out = makeFinalFigures(opts)
% Examples:
%   makeFinalFigures
%   makeFinalFigures(plots="all")
%   makeFinalFigures(3, plots=["GRdistributions","DUSP1DistPlots"])
arguments
    opts.finalModel = 2
    opts.GRmodel = 1
    opts.plots = "all"
    opts.closeAllFirst (1,1) logical = false
    opts.hideFigures = false
end

selectedPlots = normalizePlotSelection(opts.plots);
% disp("Selected plots:")
% disp(selectedPlots)
wantGRDistributions = any(strcmpi(selectedPlots, 'GRdistributions'));
wantGRMeans = any(strcmpi(selectedPlots, 'GRmeans'));
wantTotalGRMeans = any(strcmpi(selectedPlots, 'TotalGR'));
wantGRMeansAndVars = any(strcmpi(selectedPlots, 'GRmeansAndVars'));
wantTriptolide = any(strcmpi(selectedPlots, 'triptolide'));
wantAllGR = any(strcmpi(selectedPlots, 'GRdistributions'));
wantDUSP1DistPlots = any(strcmpi(selectedPlots, 'DUSP1DistPlots'));
wantGRModelComparison = any(strcmpi(selectedPlots, 'GRModelComparison'));
wantMHGR = any(strcmpi(selectedPlots, 'MHGR'));
wantMHCyt = any(strcmpi(selectedPlots, 'MHCyt'));
wantMHMech = any(strcmpi(selectedPlots, 'MHMech'));
wantGRMeansAndVarsOnly = any(strcmpi(selectedPlots,'GRModelComparison'));

if opts.hideFigures
    figVis = 'off';
    openFigVis = 'invisible';
else
    figVis = 'on';
    openFigVis = 'visible';
end

% The combined all-measurements panel reuses the GR mean/variance figures:
buildGRMeansHelper = wantGRMeans || wantTotalGRMeans || wantGRMeansAndVars || wantGRModelComparison;

out = struct();
finalModel = opts.finalModel;
GRmodel = opts.GRmodel;
out.figures = struct();

savedFullParSet = {'FinalModel_allButGR',...  %2
    'Pars_GR_DUSP1_TTP_Model_All',... %10
    'FinalModelNoSaturation',...       %3
    'FinalModelNoTimeVariation',...    %4
    'Pars_first_order_decay'};         %8

finalFullModel = 'FinalModel_allButGR';

modelChoice = normalizeFinalModelChoice(finalModel, savedFullParSet);

% ========= CENTRALIZED PALETTE (EDIT THESE) =========
C = struct();
C.NucGR_data      = [0.906, 0.529, 0.169];  % [231,135,43]/255
C.NucGR_data      = [0.906, 0.529, 0.169];  % [231,135,43]/255
C.CytoGR_data     = [0.192, 0.494, 0.760];  % approx [49,126,194]/255
C.TotalGR_data    = [0.35 0.35 0.35];
C.NucDUSP1_data   = [0.769, 0.243, 0.588];
C.CytoDUSP1_data  = [0.510, 0.235, 0.706];

C.Model_line      = [0, 0, 0];              % black for model lines
C.Shade_blend     = 0.85;                   % 0=use data color, 1=white (tint for patches)
C.TS_data         = [0.36 0.66 0.8];   % orange marker edge
% C.TS_marker_face  = [0, 0, 0];               % black marker fill
C.TS_model        = [0.122, 0.471, 0.706];   % blue model line
% C.TS_shade        = [0.80, 0.93, 0.83];      % pale mint shade

smallFont = 8;
regFont = 10;

%% Make Titration Plots
if wantDUSP1DistPlots&&~isempty(modelChoice)
    oldFigNames= {'Dusp1_SSA_Dex_Titration_TS_FracActive',...
        'Dusp1_SSA_Dex_Titration_TS_MeanActive',...
        'Dusp1_SSA_Dex_Titration_Cyt',...
        'Dusp1_SSA_Dex_Titration_Nuc'};

    f = figure;
    
    for i = 1:4
        oldFig = openfig(['savedFigures/RawFigs/',modelChoice,'/',oldFigNames{i}],'visible');
        set(oldFig,"Position",[ 1030         674         246         228]);

        ax = gca;
        grid on;
        ax.XLabel.Visible = 'off';
        ax.YLabel.Visible = 'off';
        ax.Title.Visible = 'off';
        ax.FontWeight = 'bold';
        ax.FontSize = regFont;
        switch i
            case 1
                formatLines(ax,C.TS_data);
                ax.Children(2).Color = C.TS_model;
                grid(ax,"on");
                saveas(oldFig,['savedFigures/',modelChoice,'/Dex_Titration_TS_FracActive.fig'],'fig');

                if strcmpi(modelChoice,finalFullModel)
                    saveas(oldFig,'savedFigures/FinalFigures/MainFigures/Fig4CL_Dex_Titration_TS_MeanPerActiveTS.fig','fig');
                end
                figure(f);subplot(1,3,1);
                fax = gca;
                copyobj(ax.Children,fax);
                grid(fax,"on");
                fax.XScale = 'log';
                fax.Title.String = 'TS Fraction';
                fax.XLabel.String = 'Dex Concentration (nM)';
                

            case 2
                ax = formatLines(ax,C.TS_data);
                grid(ax,"on");
                ax.Children(2).Color = C.TS_model;
                saveas(oldFig,['savedFigures/',modelChoice,'/Dex_Titration_TS_MeanPerActiveTS.fig'],'fig');
            case 3
                ax = formatLines(ax,C.CytoDUSP1_data);
                grid(ax,"on");
                saveas(oldFig,['savedFigures/',modelChoice,'/Dex_Titration_CytDUSP1.fig'],'fig');
                if strcmpi(modelChoice,finalFullModel)
                    saveas(oldFig,'savedFigures/FinalFigures/MainFigures/Fig4CR_Dex_Titration_CytDUSP1.fig','fig');
                end
                figure(f);subplot(1,3,3);
                fax = gca;
                grid(fax,"on");
                copyobj(ax.Children,fax);
                fax.XScale = 'log';
                fax.Title.String = 'DUSP1 Cyt';
                fax.XLabel.String = 'Dex Concentration (nM)';

            case 4
                ax = formatLines(ax,C.NucDUSP1_data);
                grid(ax,"on");
                saveas(oldFig,['savedFigures/',modelChoice,'/Dex_Titration_NucDUSP1.fig'],'fig');
                if strcmpi(modelChoice,finalFullModel)
                    saveas(oldFig,'savedFigures/FinalFigures/MainFigures/Fig4CM_Dex_Titration_NucDUSP1.fig','fig');
                end
                figure(f);subplot(1,3,2);
                fax = gca;
                copyobj(ax.Children,fax);
                grid(fax,"on");
                fax.XScale = 'log';
                fax.Title.String = 'DUSP1 Nuc';
                fax.XLabel.String = 'Dex Concentration (nM)';
        end
    end
    set(f,'Position',[135   865   965   273]);
    saveas(f,['savedFigures/FinalFigures/SupplementalFigures/Dex_Titration_Model_',modelChoice,'.fig'],'fig');
end
close all

%%  Distribution Plots
% Distribution plots are made one row per condition.
if opts.closeAllFirst
    close all
end

% Specify Dex concentrations
dexConc = {'100','10','1'};

% Specify the Nucleus to Cytoplasm ratio
ratio_N2C = 0.5514;

% Names of GR Models
GRModelNames = {'combinedGRModel','combinedGRModel_CytDegOnly','combinedGRModel_NucDegOnly','combinedGRModel_DexExport'};

%% Compare Cyt Models


%% GR Distributions
if wantGRDistributions
    colModel = [0,0,0];
    splitReps = false;

    for i = 1:3
        dex = dexConc{i};
        
        oldFig = openfig(['savedFigures/RawFigs/',GRModelNames{GRmodel},'_Dex',dex,'_2.fig'],openFigVis);
        for j = 1:2
            newFig = figure('Visible',figVis);
            if opts.hideFigures
                newFig.Visible = 'off';
            else
                newFig.Visible = 'on';
            end
            set(newFig,'Position',[1000        1039        1361         199]);

            if j==1
                newFig.Name = ['Nuc GR (Dex = ',dex,' nM'];
                ratio = ratio_N2C; % average ratio of nucleus to cytoplasm
                xlimv = [0,65];
                xticksv = 0:20:60;
                if splitReps
                    hidePlts = [1,2,5,6,9,10]+2;
                    showPlts = [1,2,5,6,9,10];
                    modelPlt = 1;
                else
                    hidePlts = [1,2]+2;
                    showPlts = [1,2];
                    modelPlt = 1;
                end
                dataColor = C.NucGR_data;  % nucleus data color
            else
                newFig.Name = ['Cyt GR (Dex = ',dex,' nM'];
                ratio = 1; % average ratio of cytoplasm to cytoplasm ratio
                xlimv = [0,35];
                xticksv = [0,10,20,30];
                if splitReps
                    hidePlts = [1,2,5,6,9,10];
                    showPlts = [1,2,5,6,9,10]+2;
                    modelPlt = 3;
                else
                    hidePlts = [1,2];
                    showPlts = [1,2]+2;
                    modelPlt = 3;
                end
                dataColor = C.CytoGR_data; % cytoplasm data color
            end

            subplotHandles = findobj(oldFig, 'type', 'axes');

            orderset = [7,6,5,4,3,2,1];
            for i2 = 1:length(orderset)
                isph = orderset(i2);
                subplotHandle = subplotHandles(isph);

                % Create new subplot in the new figure
                ax = subplot(1, 7, i2, 'Parent', newFig);

                % Copy contents of the original subplot to the new subplot
                copyobj(allchild(subplotHandle), ax);

                if i2~=1
                    ylabel('')
                    set(gca,'yticklabels',[]);
                end

                h = gca;

                grid on

                % Hide and show the requested lines
                for ih = hidePlts
                    if ih <= numel(h.Children)
                        h.Children(ih).Visible = 'on';
                    end
                end
                for ih = showPlts
                    if ih <= numel(h.Children)
                        h.Children(ih).Visible = 'on';
                        h.Children(ih).LineWidth = 3;
                        % Color visible data curves with the chosen data color
                        if isprop(h.Children(ih),'Color')
                            h.Children(ih).Color = dataColor;
                        end
                    end
                end

                % Model line in black
                if modelPlt <= numel(h.Children) && isprop(h.Children(modelPlt),'Color')
                    h.Children(modelPlt).Color = C.Model_line;
                    h.Children(modelPlt).LineWidth = 3;
                end

                % Adjust concentration scale for Nuclear GR.
                for ich=1:length(h.Children)
                    if isprop(h.Children(ich),'XData')
                        h.Children(ich).XData = h.Children(ich).XData/ratio;
                    end
                end

                set(gca,'xlim',xlimv,'ylim',[0,0.25]);
                set(gca,'FontSize',16,'XTick',xticksv);
            end
        end
    end
end

%% GR Means and Variances.
if buildGRMeansHelper
    clear newFig;
    for idex = 1:3
        dex = dexConc{idex};
        oldFig = openfig(['savedFigures/RawFigs/',GRModelNames{GRmodel},'_Dex',dex,'_3.fig'],openFigVis);

        for j = 1:2
            newFig{idex,j} = figure;
            if opts.hideFigures
                newFig{idex,j}.Visible = 'off';
            else
                newFig{idex,j}.Visible = 'on';
            end
            set(gcf,'Position',[1000 770   375  247]);

            copyobj(get(oldFig,'children'), newFig{idex,j});
            switch j
                case 1
                    newFig{idex,j}.Name = ['Nuc. GR vs. time (Dex = ',dex,' nM)'];
                    Plts = [1,3,4,2,5,6];
                    ratio = ratio_N2C;
                    ylimv  =[0,60];
                    dataColor = C.NucGR_data;
                case 2
                    newFig{idex,j}.Name = ['Cyt. GR vs. time (Dex = ',dex,' nM)'];
                    Plts = [2,5,6,1,3,4];
                    ratio = 1;
                    ylimv  =[0,30];
                    dataColor = C.CytoGR_data;
            end

            h = gca;
            grid on
            set(h,'Children',h.Children(Plts));
            selectedChildren = h.Children(1:min(3,numel(h.Children)));
            otherChildren = h.Children(4:end);
            set(selectedChildren, 'Visible', 'on');
            set(otherChildren, 'Visible', 'off');

            h.Children(1).LineWidth = 4;
            h.Children(2).LineWidth = 4;

            % Adjust concentration scale for Nuclear GR.
            for ich=1:3
                h.Children(ich).YData = h.Children(ich).YData/ratio;
            end
            if isprop(h.Children(1),'YNegativeDelta')
                h.Children(1).YNegativeDelta = h.Children(1).YNegativeDelta/ratio;
            end
            if isprop(h.Children(1),'YPositiveDelta')
                h.Children(1).YPositiveDelta = h.Children(1).YPositiveDelta/ratio;
            end

            legend('','','','Model \mu \pm \sigma','Model \mu','Data \mu \pm \sigma ')
            set(gca,'fontsize',16,'ylim',ylimv);
            xlabel('');
            title('');
            ylabel('');

            % Color: data errorbars = dataColor, model line = black
            formatLines(gca, dataColor);
            recolorAxes(gca, struct( ...
                'dataColor',    dataColor, ...
                'modelLine',    C.Model_line, ...
                'shadeFromModel', false));
        end
    end

    if wantGRMeans
        out.figures.GRMeans = newFig;
    end
end

%% Total GR Means vs. time
if wantTotalGRMeans

    TotalGRFig = figure('Visible',figVis);
    if opts.hideFigures
        TotalGRFig.Visible = 'off';
    else
        TotalGRFig.Visible = 'on';
    end
    set(TotalGRFig, 'Position', [1000 500 900 280]);
    TotalGRFig.Name = 'Total GR vs. time';

    rows = 1;
    cols = 3;

    leftMargin = 0.08;
    bottomMargin = 0.18;
    topMargin = 0.12;
    width = 0.94*(1-leftMargin)/cols;
    colsep = 0.04*(1-leftMargin)/cols;
    heightB = 1 - bottomMargin - topMargin;
    fontsize = 13;

    clear axsTotalGR

    for c = 1:cols
        axsTotalGR{c} = axes('Position', ...
            [(c-1)*(width+colsep)+leftMargin, bottomMargin, width, heightB]);

        fNuc = newFig{4-c,1};
        fCyt = newFig{4-c,2};

        axsTotalGR{c} = plotTotalGRFromMeanFigures( ...
            axsTotalGR{c}, fNuc, fCyt, ratio_N2C, C);

        set(axsTotalGR{c}, ...
            'XLim', [-10,190], ...
            'XTick', 0:30:180, ...
            'YLim', [14,41], ...
            'YTick', 15:5:40, ...
            'Box', 'on', ...
            'YGrid', 'on', ...
            'XGrid', 'on');

        xlabel(axsTotalGR{c}, 'Time (min)', 'FontSize', fontsize);

        if c == 1
            ylabel(axsTotalGR{c}, 'GR Total', 'FontSize', fontsize);
        else
            set(axsTotalGR{c}, 'YTickLabels', []);
        end

        title(axsTotalGR{c}, ['Dex = ', dexConc{4-c}, ' nM'], 'FontSize', fontsize);
    end

    set(TotalGRFig, 'Color', 'w');

    out.figures.TotalGR = TotalGRFig;
end

%% Compare all GR model variants at 100 nM Dex
if wantGRModelComparison
    if opts.hideFigures
        GRCompareFig = figure('Visible','off');
    else
        GRCompareFig = figure('Visible','on');
    end
    set(GRCompareFig, 'Position', [1000 350 720 850], 'Color', 'w');
    drawnow

    % Custom colors for GR model comparison plots
    clrData      = [0 0 0];             % black
    clrNucCyt    = [0.494 0.184 0.556]; % purple
    clrNucOnly   = [0.96,0.47,0.16];    % orange
    clrCytOnly   = [0.15,0.55,0.87];    % blue
    clrDexExport = [0.466 0.674 0.188]; % green

    grModelColors = { ...
        clrNucCyt, ...
        clrNucOnly, ...
        clrCytOnly, ...
        clrDexExport};

    grModelFiles = { ...
        'combinedGRModel_Dex100_3.fig', ...
        'combinedGRModel_NucDegOnly_Dex100_3.fig', ...
        'combinedGRModel_CytDegOnly_Dex100_3.fig', ...
        'combinedGRModel_DexExport_Dex100_3.fig'};

    grModelLabels = { ...
        'Nuc+Cyt Deg.', ...
        'Nuc Deg. Only', ...
        'Cyt Deg. Only', ...
        'Dex Export'};

    % Use four visually distinct line styles. Colors can be edited if
    % desired:
    grLineSpecs = {'-', '--', ':', '-.'};

    GRCompareFig = figure('Visible',figVis);
    if opts.hideFigures
        GRCompareFig.Visible = 'off';
    else
        GRCompareFig.Visible = 'on';
    end
    GRCompareFig.Name = 'GR model comparison at 100 nM Dex';
    set(GRCompareFig, 'Position', [1000 350 720 850], 'Color', 'w');

    rows = 3;
    cols = 1;

    leftMargin = 0.12;
    rightMargin = 0.04;
    bottomMargin = 0.08;
    topMargin = 0.05;
    rowsep = 0.06;

    width = 1 - leftMargin - rightMargin;
    heightB = (1 - bottomMargin - topMargin - (rows-1)*rowsep) / rows;

    axsGRCompare = cell(rows,1);

    for r = 1:rows
        axsGRCompare{r} = axes('Position', ...
            [leftMargin, ...
            bottomMargin + (rows-r)*(heightB+rowsep), ...
            width, heightB]);
        hold(axsGRCompare{r}, 'on');
        box(axsGRCompare{r}, 'on');
        grid(axsGRCompare{r}, 'on');
    end

    % ---------------------------------------------------------------------
    % Open each model's 100 nM GR mean/variance figure and extract:
    %   cyt data + model mean
    %   nuc data + model mean
    % Then compute total GR from cyt/nuc extracted values.
    % ---------------------------------------------------------------------

    for iGR = 1:numel(grModelFiles)

        figPath = fullfile('savedFigures/RawFigs', grModelFiles{iGR});

        if ~exist(figPath, 'file')
            warning('Could not find %s. Skipping this GR model.', figPath);
            continue
        end

        oldFig = openfig(figPath, openFigVis);

        [cytData, cytModel] = extractGRMeanVarFromFig(oldFig, 'cyt', ratio_N2C);
        [nucData, nucModel] = extractGRMeanVarFromFig(oldFig, 'nuc', ratio_N2C);

        % Cytoplasmic GR row.
        if iGR == 1
            errorbar(axsGRCompare{1}, ...
                cytData.x, cytData.y, cytData.yneg, cytData.ypos, ...
                'o', ...
                'Color', clrData, ...
                'MarkerFaceColor', clrData, ...
                'MarkerEdgeColor', clrData, ...
                'LineWidth', 1.4, ...
                'CapSize', 8, ...
                'DisplayName', 'Data');
        end

        plot(axsGRCompare{1}, cytModel.x, cytModel.y, ...
            grLineSpecs{iGR}, ...
            'Color', grModelColors{iGR}, ...
            'LineWidth', 2.5, ...
            'DisplayName', grModelLabels{iGR});

        % Nuclear GR row.
        if iGR == 1
            errorbar(axsGRCompare{2}, ...
                nucData.x, nucData.y, nucData.yneg, nucData.ypos, ...
                'o', ...
                'Color', clrData, ...
                'MarkerFaceColor', clrData, ...
                'MarkerEdgeColor', clrData, ...
                'LineWidth', 1.4, ...
                'CapSize', 8, ...
                'DisplayName', 'Data');
        end

        plot(axsGRCompare{2}, nucModel.x, nucModel.y, ...
            grLineSpecs{iGR}, ...
            'Color', grModelColors{iGR}, ...
            'LineWidth', 2.5, ...
            'DisplayName', grModelLabels{iGR});

        % Total GR row.
        totalData.x = cytData.x;
        totalData.y = cytData.y + ratio_N2C*nucData.y;
        totalData.yneg = sqrt(cytData.yneg.^2 + (ratio_N2C*nucData.yneg).^2);
        totalData.ypos = sqrt(cytData.ypos.^2 + (ratio_N2C*nucData.ypos).^2);

        totalModel.x = cytModel.x;
        totalModel.y = cytModel.y + ratio_N2C*nucModel.y;

        if iGR == 1
            errorbar(axsGRCompare{3}, ...
                totalData.x, totalData.y, totalData.yneg, totalData.ypos, ...
                'o', ...
                'Color', clrData, ...
                'MarkerFaceColor', clrData, ...
                'MarkerEdgeColor', clrData, ...
                'LineWidth', 1.4, ...
                'CapSize', 8, ...
                'DisplayName', 'Data');
        end

        plot(axsGRCompare{3}, totalModel.x, totalModel.y, ...
            grLineSpecs{iGR}, ...
            'Color', grModelColors{iGR}, ...
            'LineWidth', 2.5, ...
            'DisplayName', grModelLabels{iGR});

        close(oldFig);
    end

    ylabel(axsGRCompare{1}, 'Cyt. GR', 'FontSize', 14);
    ylabel(axsGRCompare{2}, 'Nuc. GR', 'FontSize', 14);
    ylabel(axsGRCompare{3}, 'Total GR', 'FontSize', 14);
    xlabel(axsGRCompare{3}, 'Time (min)', 'FontSize', 14);

    title(axsGRCompare{1}, 'Cytoplasmic GR, Dex = 100 nM', 'FontSize', 14);
    title(axsGRCompare{2}, 'Nuclear GR, Dex = 100 nM', 'FontSize', 14);
    title(axsGRCompare{3}, 'Total GR, Dex = 100 nM', 'FontSize', 14);

    set(axsGRCompare{1}, 'XLim', [-10 190], 'YLim', [0 30], ...
        'XTick', 0:30:180, 'FontSize', 13);
    set(axsGRCompare{2}, 'XLim', [-10 190], 'YLim', [0 62], ...
        'XTick', 0:30:180, 'FontSize', 13);
    set(axsGRCompare{3}, 'XLim', [-10 190], 'YLim', [14 41], ...
        'XTick', 0:30:180, 'FontSize', 13);

    set(axsGRCompare{1}, 'XTickLabels', []);
    set(axsGRCompare{2}, 'XTickLabels', []);

    legend(axsGRCompare{1}, 'Location', 'best', 'FontSize', 11);

    out.figures.GRModelComparison = GRCompareFig;

    if ~exist('savedFigures', 'dir')
        mkdir('savedFigures');
    end

    savefig(GRCompareFig, fullfile('savedFigures/GRModels', 'GRModelComparison_Dex100.fig'));
    saveas(GRCompareFig, fullfile('savedFigures/GRModels', 'GRModelComparison_Dex100.png'));
    if GRmodel==1
        set(GRCompareFig,'Position',[ 374   718   200   467])
        for i = 1:length(GRCompareFig.Children)
            % GRCompareFig.Children(i).FontSize = 10;
            if ~strcmpi(GRCompareFig.Children(i).Type,'legend')
                GRCompareFig.Children(i).FontSize = regFont;
                GRCompareFig.Children(i).FontWeight = 'bold';
                GRCompareFig.Children(i).XLabel.Visible = 'off';
                GRCompareFig.Children(i).YLabel.Visible = 'off';
            else
                GRCompareFig.Children(i).Visible = 'off';
            end
        end
        GRCompareFig.Children(1).XTickLabelRotation = 0;
        saveas(GRCompareFig,'savedFigures/FinalFigures/MainFigures/Fig2C_GRModels_Comparison.fig','fig');
    end
end

%% GR Model Fits Means and Variances
if wantGRMeansAndVarsOnly

    AllGRDataFig = figure('Visible',figVis);
    set(AllGRDataFig,'Position',[1200         618         725         450]);

    % Define the number of rows and columns
    rows=4;
    cols=3;

    leftMargin = 0.08;
    bottomMargin = 0.07;
    topMargin = 0.02;
    width = 0.94*(1-leftMargin)/cols;
    colsep = 0.04*(1-leftMargin)/cols;
    heightB = 0.91*(1-bottomMargin-topMargin)/rows;
    rowsep = 0.08*(1-bottomMargin)/rows;
    fontsize = 13;

    rowLabels = {'Dex (nM)','GR Cyt (AU)','GR Nuc (AU)','GR Total'};
    % Create axes in a loop
    for r = 1:rows
        for c = 1:cols
            % Calculate position for each axes
            axs{r,c} = axes('Position', [(c-1)*(width+colsep)+leftMargin, (rows-r)*(rowsep+heightB)+bottomMargin, width, heightB]);

            if r==1
                plot([-20,0,0,200],(log10(eval(dexConc{4-c}))+1)*[0,0,1,1]+1,'k','LineWidth',2);
                ylimv = [0,5];
                yticks = 0:1:4;
                set(gca,'YTickLabel',{'','0','1','10','100',''})
                if r==1
                    ti = title(['Dex = ',dexConc{4-c},' nM']);
                    ti.FontSize = fontsize;
                end

            elseif r==2 % Cyto GR
                ax1 = newFig{4-c,2}.Children(2);
                children = get(ax1, 'Children'); % Get all children of ax1
                copyobj(children, axs{r,c}); % Copy all children to ax2
                axs{r,c} = formatLines(axs{r,c}, C.CytoGR_data);   % color data errorbars + standardize
                recolorAxes(axs{r,c}, struct('dataColor', C.CytoGR_data, 'modelLine', C.Model_line, 'shadeFromModel', false));
                ylimv = [0,30];
                yticks = 0:10:30;

            elseif r==3 % Nuc GR
                ax1 = newFig{4-c,1}.Children(2);
                children = get(ax1, 'Children'); % Get all children of ax1
                copyobj(children, axs{r,c}); % Copy all children to ax2
                axs{r,c} = formatLines(axs{r,c}, C.NucGR_data);
                recolorAxes(axs{r,c}, struct('dataColor', C.NucGR_data, 'modelLine', C.Model_line, 'shadeFromModel', false));
                ylimv = [0,62];
                yticks = 0:15:60;

            elseif r==4 % Total GR = cytoplasmic + ratio*nuclear

                fNuc = newFig{4-c,1};
                fCyt = newFig{4-c,2};

                axs{r,c} = plotTotalGRFromMeanFigures( ...
                    axs{r,c}, fNuc, fCyt, ratio_N2C, C);

                ylimv = [14,41];
                yticks = 15:5:40;
          
            else
                ylimv = [0,1];
                yticks = 0:0.2:1;
            end

            % Set x-ticks and y-ticks
            set(axs{r,c}, 'XLim',[-10,190],'XTick', 0:30:180);
            if r == rows % Only bottom row gets x-ticks
                xlabel('Time (min)','FontSize',fontsize);
            else
                set(axs{r,c}, 'XTickLabels', []); % Remove x-ticks for other rows
            end

            set(axs{r,c}, 'YTick', yticks,'YLim',ylimv);
            if c == 1 % Only leftmost column gets y-ticks
                ylab = ylabel(rowLabels{r},'FontSize',13);
                ylab.Position(1) = -35;
            else
                set(axs{r,c}, 'YTickLabels', []); % Remove y-ticks for other columns
            end

            % Minimize space between axes
            set(axs{r,c}, 'Box', 'on', 'YGrid', 'on', 'XGrid', 'on'); % Optional: add box around axes
        end
    end

    % Adjust figure properties to minimize space
    set(gcf, 'Color', 'w'); % Set figure background color

    saveas(AllGRDataFig,['savedFigures/GRModels/',GRModelNames{GRmodel},'_GRvsTime.fig'],'fig');


    out.figures.AllGRDataFig = AllGRDataFig;

    AllGRDataFig = openfig(['savedFigures/GRModels/',GRModelNames{GRmodel},'_GRvsTime.fig']);
    set(AllGRDataFig,'Position',[323   753   497   487]);

    for i = 1:length(AllGRDataFig.Children)
        AllGRDataFig.Children(i).FontSize = regFont;
        AllGRDataFig.Children(i).FontWeight ='bold';
        if GRmodel==1
            AllGRDataFig.Children(i).XLabel.Visible = 'off';
            AllGRDataFig.Children(i).YLabel.Visible = 'off';
        end
    end
    AllGRDataFig.Children(1).Color = [0 0 0]*0.15+.85;
    AllGRDataFig.Children(2).Color = [0 0 0]*0.1+.9;
    AllGRDataFig.Children(3).Color = [0 0 0]*0.05+.95;
    AllGRDataFig.Children(3).YTick = [15,25,35];
    AllGRDataFig.Children(4).Color = [1 0 0]*0.15+.85;
    AllGRDataFig.Children(5).Color = [1 0 0]*0.1+.9;
    AllGRDataFig.Children(6).Color = [1 0 0]*0.05+.95;
    AllGRDataFig.Children(7).Color = [0 0 1]*0.15+.85;
    AllGRDataFig.Children(8).Color = [0 0 1]*0.1+.9;
    AllGRDataFig.Children(9).Color = [0 0 1]*0.05+.95;
    hold(AllGRDataFig.Children(10),'on');
    fill(AllGRDataFig.Children(10),[0,200,200,0],[1,1,4,4],[0.85 0.85 0.85]);
    hold(AllGRDataFig.Children(11),'on');
    fill(AllGRDataFig.Children(11),[0,200,200,0],[1,1,3,3],[0.9 0.9 0.9]);
    hold(AllGRDataFig.Children(12),'on');
    fill(AllGRDataFig.Children(12),[0,200,200,0],[1,1,2,2],[0.95 0.95 0.95]);

    if GRmodel==1
        saveas(AllGRDataFig,'savedFigures/FinalFigures/MainFigures/Fig2E_GR_MeansAndVars_v_Dex_and_Time.fig','fig');
    else

        for i = 1:3
            AllGRDataFig.Children(i).XTickLabelRotation = 0;
            AllGRDataFig.Children(i).XRuler.TickLabelGapOffset = -1;    % distance from axis (R2023a+)
        end
        AllGRDataFig.Children(1).XLabel.Visible = 'off';
        AllGRDataFig.Children(3).XLabel.Visible = 'off';
        AllGRDataFig.Children(2).XLabel.Position(2) = 9;
        AllGRDataFig.Children(2).XLabel.String = '- - - - - - Time (min) - - - - - -';
        saveas(AllGRDataFig,['savedFigures/FinalFigures/SupplementalFigures/GR_MeansAndVars_AltModel_',GRModelNames{GRmodel},'.fig'],'fig');
    end
end


%% All GR and DUSP1 means and variances vs time Dex = 1,10,100 nm
if wantGRMeansAndVars&&~isempty(modelChoice)&&GRmodel==1
    % Create a figure
    version = '';

    AllDataFig = figure('Visible',figVis);
    set(AllDataFig,'Position',[1200         618         725         900]);

    % Define the number of rows and columns  
    rows=8;
    cols = 3;

    leftMargin = 0.08;
    bottomMargin = 0.05;
    topMargin = 0.02;
    width = 0.94*(1-leftMargin)/cols;
    colsep = 0.04*(1-leftMargin)/cols;
    heightB = 0.91*(1-bottomMargin-topMargin)/rows;
    rowsep = 0.08*(1-bottomMargin)/rows;
    fontsize = 13;

    rowLabels = {'Dex (nM)','GR Cyt (AU)','GR Nuc (AU)','GR Total','TS Fraction','RNA per TS','DUSP1-Nuc','DUSP1-Cyt'};
    % Create axes in a loop
    for r = 1:rows
        for c = 1:cols
            % Calculate position for each axes
            axs{r,c} = axes('Position', [(c-1)*(width+colsep)+leftMargin, (rows-r)*(rowsep+heightB)+bottomMargin, width, heightB]);

            if r==1
                plot([-20,0,0,200],(log10(eval(dexConc{4-c}))+1)*[0,0,1,1]+1,'k','LineWidth',2);
                ylimv = [0,5];
                yticks = 0:1:4;
                set(gca,'YTickLabel',{'','0','1','10','100',''});
                if r==1
                    ti = title(['Dex = ',dexConc{4-c},' nM']);
                    ti.FontSize = fontsize;
                end

            elseif r==2 % Cyto GR
                ax1 = newFig{4-c,2}.Children(2);
                children = get(ax1, 'Children'); % Get all children of ax1
                copyobj(children, axs{r,c}); % Copy all children to ax2
                axs{r,c} = formatLines(axs{r,c}, C.CytoGR_data);   % color data errorbars + standardize
                recolorAxes(axs{r,c}, struct('dataColor', C.CytoGR_data, 'modelLine', C.Model_line, 'shadeFromModel', false));
                ylimv = [0,30];
                yticks = 0:10:30;

            elseif r==3 % Nuc GR
                ax1 = newFig{4-c,1}.Children(2);
                children = get(ax1, 'Children'); % Get all children of ax1
                copyobj(children, axs{r,c}); % Copy all children to ax2
                axs{r,c} = formatLines(axs{r,c}, C.NucGR_data);
                recolorAxes(axs{r,c}, struct('dataColor', C.NucGR_data, 'modelLine', C.Model_line, 'shadeFromModel', false));
                ylimv = [0,62];
                yticks = 0:15:60;

            elseif r==4 % Total GR = cytoplasmic + ratio*nuclear

                fNuc = newFig{4-c,1};
                fCyt = newFig{4-c,2};

                axs{r,c} = plotTotalGRFromMeanFigures( ...
                    axs{r,c}, fNuc, fCyt, ratio_N2C, C);

                ylimv = [14,41];
                yticks = 15:5:40;

            elseif r==5 % TS fraction
                f = open(['savedFigures/RawFigs/',modelChoice,'/Dusp1_TS_',dexConc{4-c},'_MeansAndFracs.fig']);
                ax1 = f.Children(1);
                children = get(ax1, 'Children'); % Get all children of ax1
                children(3).FaceColor = [0.9,1,0.9];
                % figure(AllDataFig);
                set(groot,'CurrentFigure',AllDataFig);
                copyobj(children, axs{r,c}); % Copy all children to ax2
                ylimv = [0,0.55];
                yticks = 0:0.1:0.5;
                axs{r,c} = formatLines(axs{r,c}, C.TS_data); % no recolor; just formatting

            elseif r==6 % RNA per TS
                f = open(['savedFigures/RawFigs/',modelChoice,'/Dusp1_TS_',dexConc{4-c},'_MeansAndFracs.fig']);
                ax1 = f.Children(2);
                children = get(ax1, 'Children'); % Get all children of ax1
                children(3).FaceColor = [0.9,1,0.9];
                set(groot,'CurrentFigure',AllDataFig);
                % figure(AllDataFig);
                copyobj(children, axs{r,c}); % Copy all children to ax2
                ylimv = [0,20];
                yticks = 0:5:20;
                axs{r,c} = formatLines(axs{r,c}, C.TS_data); % no recolor; just formatting

            elseif r==7 % DUSP1 Nuc (SSA Mn/Var)
                f = open(['savedFigures/RawFigs/',modelChoice,'/Dusp1_SSA_',dexConc{4-c},'_NucMnVar.fig']);
                ax1 = f.Children;
                children = get(ax1, 'Children'); % Get all children of ax1
                set(groot,'CurrentFigure',AllDataFig);
                % figure(AllDataFig);
                copyobj(children, axs{r,c}); % Copy all children to ax2
                ylimv = [0,110];
                yticks = 0:20:110;
                axs{r,c} = formatLines(axs{r,c}, C.NucDUSP1_data);
                recolorAxes(axs{r,c}, struct('dataColor', C.NucDUSP1_data, 'modelLine', C.Model_line, 'shadeFromModel', false));

            elseif r==8 % DUSP1 Cyt (SSA Mn/Var)
                f = open(['savedFigures/RawFigs/',modelChoice,'/Dusp1_SSA_',dexConc{4-c},'_CytMnVar.fig']);
                ax1 = f.Children;
                children = get(ax1, 'Children'); % Get all children of ax1
                set(groot,'CurrentFigure',AllDataFig);
                % figure(AllDataFig);
                copyobj(children, axs{r,c}); % Copy all children to ax2
                ylimv = [0,210];
                yticks = 0:50:200;
                axs{r,c} = formatLines(axs{r,c}, C.CytoDUSP1_data);
                recolorAxes(axs{r,c}, struct('dataColor', C.CytoDUSP1_data, 'modelLine', C.Model_line, 'shadeFromModel', false));

            else
                ylimv = [0,1];
                yticks = 0:0.2:1;
            end

            % Set x-ticks and y-ticks
            set(axs{r,c}, 'XLim',[-10,190],'XTick', 0:30:180);
            if r == rows % Only bottom row gets x-ticks
                xlabel('Time (min)','FontSize',fontsize);
            else
                set(axs{r,c}, 'XTickLabels', []); % Remove x-ticks for other rows
            end

            set(axs{r,c}, 'YTick', yticks,'YLim',ylimv);
            if c == 1 % Only leftmost column gets y-ticks
                ylab = ylabel(rowLabels{r},'FontSize',13);
                ylab.Position(1) = -35;
            else
                set(axs{r,c}, 'YTickLabels', []); % Remove y-ticks for other columns
            end

            % Minimize space between axes
            set(axs{r,c}, 'Box', 'on', 'YGrid', 'on', 'XGrid', 'on'); % Optional: add box around axes
        end
    end

    % Adjust figure properties to minimize space
    set(gcf, 'Color', 'w'); % Set figure background color
    saveas(AllDataFig,['savedFigures/',modelChoice,'/AllMsmtsVsTime.fig'],'fig');
    out.figures.AllDataFig = AllDataFig;

    AllDataFig = openfig(['savedFigures/',modelChoice,'/AllMsmtsVsTime.fig']);
    % set(AllGRDataFig,'Position',[323   753   497   487]);

    for i = 1:length(AllDataFig.Children)
        AllDataFig.Children(i).FontSize = regFont;
        AllDataFig.Children(i).FontWeight ='bold';
        AllDataFig.Children(i).XLabel.Visible = 'off';
        AllDataFig.Children(i).YLabel.Visible = 'off';
    end
    AllDataFig.Children(1+12).Color = [0 0 0]*0.15+.85;
    AllDataFig.Children(2+12).Color = [0 0 0]*0.1+.9;
    AllDataFig.Children(3+12).Color = [0 0 0]*0.05+.95;
    AllDataFig.Children(3+12).YTick = [15,25,35];
    AllDataFig.Children(4+12).Color = [1 0 0]*0.15+.85;
    AllDataFig.Children(5+12).Color = [1 0 0]*0.1+.9;
    AllDataFig.Children(6+12).Color = [1 0 0]*0.05+.95;
    AllDataFig.Children(7+12).Color = [0 0 1]*0.15+.85;
    AllDataFig.Children(8+12).Color = [0 0 1]*0.1+.9;
    AllDataFig.Children(9+12).Color = [0 0 1]*0.05+.95;
    hold(AllDataFig.Children(10+12),'on');
    fill(AllDataFig.Children(10+12),[0,200,200,0],[1,1,4,4],[0.85 0.85 0.85]);
    AllDataFig.Children(10+12).Children(end)
    hold(AllDataFig.Children(11+12),'on');
    fill(AllDataFig.Children(11+12),[0,200,200,0],[1,1,3,3],[0.9 0.9 0.9]);
    hold(AllDataFig.Children(12+12),'on');
    fill(AllDataFig.Children(12+12),[0,200,200,0],[1,1,2,2],[0.95 0.95 0.95]);

    AllDataFig.Children(1+9).Color = [0 1 0]*0.15+.85;
    AllDataFig.Children(2+9).Color = [0 1 0]*0.1+.9;
    AllDataFig.Children(3+9).Color = [0 1 0]*0.05+.95;
    AllDataFig.Children(1+6).Color = [0 1 0]*0.15+.85;
    AllDataFig.Children(2+6).Color = [0 1 0]*0.1+.9;
    AllDataFig.Children(3+6).Color = [0 1 0]*0.05+.95;

    AllDataFig.Children(1+3).Color = C.NucDUSP1_data*0.15+.85;
    AllDataFig.Children(2+3).Color = C.NucDUSP1_data*0.1+.9;
    AllDataFig.Children(3+3).Color = C.NucDUSP1_data*0.05+.95;
    AllDataFig.Children(1).Color = C.CytoDUSP1_data*0.15+.85;
    AllDataFig.Children(2).Color = C.CytoDUSP1_data*0.1+.9;
    AllDataFig.Children(3).Color = C.CytoDUSP1_data*0.05+.95;

    if strcmp(modelChoice,finalFullModel)
        saveas(AllDataFig,'savedFigures/FinalFigures/MainFigures/Fig4A_GR_and_DUSP1_MeansAndVars_v_Dex_and_Time.fig','fig');
    else
        saveas(AllDataFig,['savedFigures/FinalFigures/SupplementalFigures/GR_and_DUSP1_MeansAndVars_Model_',modelChoice,'.fig'],'fig');
    end
end


%% All Means And Variances vs. time - TRiptolide
if wantTriptolide&&~isempty(modelChoice)
    % Create a figure
    TplDataFig = figure('Visible',figVis);
    set(TplDataFig,'Position',[200         618         925         700]);

    % Define the number of rows and columns
    rows = 4;
    cols = 5;

    leftMargin = 0.08;
    bottomMargin = 0.05;
    topMargin = 0.02;
    width = 0.94*(1-leftMargin)/cols;
    colsep = 0.04*(1-leftMargin)/cols;
    heightB = 0.91*(1-bottomMargin-topMargin)/rows;
    rowsep = 0.08*(1-bottomMargin)/rows;
    fontsize = 13;

    tplA = {'Tpl_0min','Tpl_20min','Tpl_75min','Tpl_150min','Tpl_180min'};
    tTrp = [0,20,75,150,180];
    tmax = [60,80,140,220,250];
    rowLabels = {'Tpl (nM)','TS Fraction','DUSP1-Nuc','DUSP1-Cyt'};

    % Create axes in a loop
    for r = 1:rows
        for c = 1:cols
            % Calculate position for each axes
            axs{r,c} = axes('Position', [(c-1)*(width+colsep)+leftMargin, (rows-r)*(rowsep+heightB)+bottomMargin, width, heightB]);

            if r==1
                plot([-20,tTrp(c),tTrp(c),240],[0,0,1,1],'k','LineWidth',2);
                ylimv = [-0.5,2];
                yticks = 0:1:2;
                ti = title(['t_{Trp} = ',num2str(tTrp(c)),' min']);
                ti.FontSize = fontsize;

            elseif r==2 % TS fraction
                f = open(['savedFigures/RawFigs/',modelChoice,'/Dusp1_TS_',tplA{c},'_MeansAndFracs.fig']);
                ax1 = f.Children(1);
                children = get(ax1, 'Children'); % Get all children of ax1
                children(3).FaceColor = [0.9,1,0.9];
                set(groot,'CurrentFigure',TplDataFig);
                copyobj(children, axs{r,c}); % Copy all children to ax2
                ylimv = [0,0.55];
                yticks = 0:0.1:0.5;
                axs{r,c} = formatLines(axs{r,c}, C.TS_data); % standardize only
                set(gca,'xlim',[0,tmax(c)]);

                % elseif r==3 % RNA per TS (commented out)
            elseif r ==3 % DUSP1 Nuc
                f = open(['savedFigures/RawFigs/',modelChoice,'/Dusp1_SSA_',tplA{c},'_NucMnVar.fig']);
                ax1 = f.Children;
                children = get(ax1, 'Children'); % Get all children of ax1
                set(groot,'CurrentFigure',TplDataFig);
                copyobj(children, axs{r,c}); % Copy all children to ax2
                ylimv = [0,110];
                yticks = 0:20:110;
                axs{r,c} = formatLines(axs{r,c}, C.NucDUSP1_data); % tint patches from data color
                recolorAxes(axs{r,c}, struct('dataColor', C.NucDUSP1_data, 'modelLine', C.Model_line, 'shadeFromModel', false));
                set(gca,'xlim',[0,tmax(c)]);

            elseif r ==4 % DUSP1 Cyt
                f = open(['savedFigures/RawFigs/',modelChoice,'/Dusp1_SSA_',tplA{c},'_CytMnVar.fig']);
                ax1 = f.Children;
                children = get(ax1, 'Children'); % Get all children of ax1
                set(groot,'CurrentFigure',TplDataFig);
                copyobj(children, axs{r,c}); % Copy all children to ax2
                ylimv = [0,210];
                yticks = 0:50:200;
                axs{r,c} = formatLines(axs{r,c}, C.CytoDUSP1_data);
                recolorAxes(axs{r,c}, struct('dataColor', C.CytoDUSP1_data, 'modelLine', C.Model_line, 'shadeFromModel', false));
                set(gca,'xlim',[0,tmax(c)]);

            else
                ylimv = [0,1];
                yticks = 0:0.2:1;
            end

            set(gca,'xlim',[-5,tmax(c)+5]);

            % Set x-ticks and y-ticks
            if r == rows % Only bottom row gets x-ticks
                xlabel('Time (min)','FontSize',fontsize);
            else
                set(axs{r,c}, 'XTickLabels', []); % Remove x-ticks for other rows
            end

            set(axs{r,c}, 'YTick', yticks,'YLim',ylimv);
            if c == 1 % Only leftmost column gets y-ticks
                ylab = ylabel(rowLabels{r},'FontSize',13);
                ylab.Position(1) = -20;
            else
                set(axs{r,c}, 'YTickLabels', []); % Remove y-ticks for other columns
            end

            % Minimize space between axes
            set(axs{r,c}, 'Box', 'on', 'YGrid', 'on', 'XGrid', 'on'); % Optional: add box around axes

            hold on;
            plot(tTrp(c)*[1,1],ylimv,'c--','LineWidth',2);
        end
    end

    % Adjust figure properties to minimize space
    set(gcf, 'Color', 'w'); % Set figure background color
    saveas(TplDataFig,['savedFigures/',modelChoice,'/TplMsmtsVsTime.fig'],'fig');

    % Adjust for saving of final figures
    f = openfig(['savedFigures/',modelChoice,'/TplMsmtsVsTime.fig']);
    axOld = findall(f,'Type','axes');
    axOldSort = axOld(end:-1:1);
    vv = [2,3];
    vvv = [vv,5+vv,10+vv,15+vv];
    axOldSubset = axOldSort(vvv);
    newFig = figure;
    set(newFig,'Position',[373   319   453   749]);
    copyobj(axOldSubset,newFig);
    for i = 1:2
        set(newFig.Children(i),'Position',axOldSort(i).Position);
        set(newFig.Children(i+2),'Position',axOldSort(i+5).Position);
        set(newFig.Children(i+4),'Position',axOldSort(i+10).Position);
        set(newFig.Children(i+6),'Position',axOldSort(i+15).Position); 
    end
    for i = 1:length(newFig.Children)
        newFig.Children(i).Position(1,3) = newFig.Children(i).Position(1,3)*5/2;
        newFig.Children(i).Position(1,1) = newFig.Children(1).Position(1,1) + (newFig.Children(i).Position(1,1)-newFig.Children(1).Position(1,1))*5/2;
        newFig.Children(i).FontSize = smallFont;
        newFig.Children(i).FontWeight = 'bold';
        newFig.Children(i).XLabel.Visible = 'off';
        newFig.Children(i).YLabel.Visible = 'off';
    end
    for i = 2:length(newFig.Children)/2
        newFig.Children((i-1)*2+1).YTickLabel = axs{i,1}.YTickLabel;
    end
    newFig.Children(1).YTick = [0,1];
    newFig.Children(1).YTickLabel = {'0','5'};

    for i = 2:2:length(newFig.Children)
        newFig.Children(i).XTick = [0:30:120];
        newFig.Children(i-1).XTick = [0:30:120];
    end

    saveas(f,['savedFigures/FinalFigures/SupplementalFigures/DUSP1_MeansAndVars_TPL_',modelChoice,'.fig'],'fig');

    if strcmp(modelChoice,finalFullModel)
        saveas(newFig,'savedFigures/FinalFigures/MainFigures/Fig4B_DUSP1_v_TPL.fig','fig');
    end

    out.figures.TplDataFig = TplDataFig;

    close all

    figfileDusp1Tpl = 'savedFigures/FinalModel_allButGR/TplMsmtsVsTime.fig';
    TPL_degradation_Plots;
    F{1} = figure(112);
    F{2} = figure(113);
    F{3} = figure(122);
    F{4} = figure(123);

    f = figure;
    for i1 = 1:4
        subplot(2,2,i1);
        ax = gca;
        copyobj(F{i1}.Children(1).Children, ax);
        ax.Children(2).Color = [0,0,0];
        switch i1
            case 1
                ax = formatLines(ax, C.NucDUSP1_data);
                ax.YLabel.String = 'Nuc Dusp1';
            case 2
                ax = formatLines(ax, C.NucDUSP1_data);
            case 3
                ax = formatLines(ax, C.CytoDUSP1_data);
                ax.YLabel.String = 'Cyt Dusp1';
                ax.XLabel.String = 'Time (min)';
            case 4
                ax = formatLines(ax, C.CytoDUSP1_data);
                ax.XLabel.String = 'Time (min)';
        end
        ax.YLim = [7e-2,3];
        ax.YTick = [0.1 0.2 0.4 0.8 1.6];
        ax.XTick = 0:10:60;
        ax.FontSize = regFont;
        ax.YScale = 'log';
        grid(ax,'minor');
        grid(ax,'on');
        ax.Title.String = ['DUSP1 Decay for Model ',modelChoice];
    end

    for i = 1:4
        set(F{i},'Position',[941   786   263   218]);
    end
    for i = [1,3]
        F{i}.Children(1).Children(2).Color = [0,0,0];
        F{i+1}.Children(1).Children(2).Color = [0,0,0];
    end
    for i = 1:2
        F{i}.Children(1) = formatLines(F{i}.Children(1), C.NucDUSP1_data);
        F{i+2}.Children(1) = formatLines(F{i+2}.Children(1), C.CytoDUSP1_data);
    end
    for i = 1:4
        F{i}.Children(1).YLim = [7e-2,3];
        F{i}.Children(1).YTick = [0.1 0.2 0.4 0.8 1.6];
        F{i}.Children(1).XTick = 0:10:60;
        F{i}.Children(1).FontSize = regFont;
        grid(F{i}.Children(1),'minor');
        grid on
    end
    if strcmp(modelChoice,finalFullModel)
        saveas(F{1},'savedFigures/FinalFigures/MainFigures/Fig4DTL_DUSP1_v_TPL_NucEarly.fig','fig');
        saveas(F{2},'savedFigures/FinalFigures/MainFigures/Fig4DTR_DUSP1_v_TP_NucLate.fig','fig');
        saveas(F{3},'savedFigures/FinalFigures/MainFigures/Fig4DBL_DUSP1_v_TPL_CytEarly.fig','fig');
        saveas(F{4},'savedFigures/FinalFigures/MainFigures/Fig4DBR_DUSP1_v_TPL_CytLate.fig','fig');
    end
    saveas(f,['savedFigures/FinalFigures/SupplementalFigures/DecayRatesModel_',modelChoice,'.fig'],'fig');
end

%% All GR distributions
if wantAllGR
    % Create a figure
    try close(AllDistsFig); catch; end
    version = '';

    AllDistsFig = figure('Visible',figVis);
    set(AllDistsFig,'Position',[566         833        1359         522]);

    % Define the number of rows and columns
    rows = 6;
    cols = 7;

    leftMargin = 0.08;
    bottomMargin = 0.05;
    topMargin = 0.02;
    width = 0.86*(1-leftMargin)/cols;
    colsep = 0.08*(1-leftMargin)/cols;
    heightB = 0.8*(1-bottomMargin-topMargin)/rows;
    rowsep = 0.18*(1-bottomMargin)/rows;
    fontsize = 13;

    dexConc = [1,10,100];
    tGRArr = [0,10,30,50,75,120,180];

    rowLabels = {{'1nM Dex';'Cytoplasm'},{'1nM Dex';'Nucleus'},{'10nM Dex';'Cytoplasm'},{'10nM Dex';'Nucleus'},{'100nM Dex';'Cytoplasm'},{'100nM Dex';'Nucleus'}};
    % Create axes in a loop
    for r = 1:rows
        oldFig = openfig(['savedFigures/RawFigs/',GRModelNames{GRmodel},'_Dex',num2str(dexConc(floor((r-1)/2)+1)),'_2.fig'],openFigVis);
        set(groot,'CurrentFigure',AllDistsFig);

        for c = 1:cols

            % Calculate position for each axes
            axs{r,c} = axes('Position', [(c-1)*(width+colsep)+leftMargin, (rows-r)*(rowsep+heightB)+bottomMargin, width, heightB]);

            if mod(r,2)==0
                copyobj(oldFig.Children(9-c).Children(1:2), axs{r,c});
                set(gca,'YLim',[0,0.2],'XLim',[-2,65]);
                ratio = ratio_N2C; %Average ratio of nucleus to cytoplasm.
                for ich = 1:2
                    axs{r,c}.Children(ich).XData = axs{r,c}.Children(ich).XData/ratio;
                end
                axs{r,c}.Children(2).Color = C.NucGR_data;
                axs{r,c}.Children(1).Color = C.Model_line;

            else
                copyobj(oldFig.Children(9-c).Children(3:4), axs{r,c});
                set(gca,'YLim',[0,0.3],'XLim',[-2,35]);
                ratio = 1.0; %Average ratio of nucleus to cytoplasm.
                axs{r,c}.Children(2).Color = C.CytoGR_data;
                axs{r,c}.Children(1).Color = C.Model_line;
            end

            if c~=1
                set(gca,'yticklabels',[]);
            else
                ylab = ylabel(rowLabels{r});
                xlim = get(gca,'XLim');
                ylab.Position(1) = -0.22*(xlim(2)-xlim(1));
            end

            if r==1
                title(['t = ',num2str(tGRArr(c)),' min']);
            end

            % % Minimize space between axes
            set(axs{r,c}, 'Box', 'on', 'YGrid', 'on', 'XGrid', 'on'); % Optional: add box around axes
        end
    end

    % Adjust figure properties to minimize space
    set(gcf, 'Color', 'w'); % Set figure background color
    % saveas(AllDataFig,'savedFigures/AllMsmtsVsTime.fig','fig');

    out.figures.AllGRDistsFig = AllDistsFig;

    saveas(AllDistsFig,['savedFigures/GRModels/',GRModelNames{GRmodel},'_Distributions.fig'],'fig');

    set(AllDistsFig,'Position',[270   703   712   522])
    for i = 1:length(AllDistsFig.Children)
        if GRmodel==1
            AllDistsFig.Children(i).XLabel.Visible = 'off';
            AllDistsFig.Children(i).YLabel.Visible = 'off';
        end
        AllDistsFig.Children(i).FontSize = smallFont;
        AllDistsFig.Children(i).FontWeight = 'bold';
        if i<=7
            AllDistsFig.Children(i).Color = (1-C.CytoGR_data)*0.15+.85;
        elseif i<=14
            AllDistsFig.Children(i).Color = (1-C.NucGR_data)*0.15+.85;
        elseif i<=21
            AllDistsFig.Children(i).Color = (1-C.CytoGR_data)*0.1+.9;
        elseif i<=28
            AllDistsFig.Children(i).Color = (1-C.NucGR_data)*0.1+.9;
        elseif i<=35
            AllDistsFig.Children(i).Color = (1-C.CytoGR_data)*0.05+.95;
        else
            AllDistsFig.Children(i).Color = (1-C.NucGR_data)*0.05+.95;
        end
    end
    if GRmodel==1
        saveas(AllDistsFig,'savedFigures/FinalFigures/MainFigures/Fig2D_GRDistributions.fig','fig');
    else
        AllDistsFig.Children(4).XLabel.String = '- - - - - - Molecular Count (AU) - - - - - -';
        saveas(AllDistsFig,['savedFigures/FinalFigures/SupplementalFigures/GR_Distributions_AltModel_',GRModelNames{GRmodel},'.fig'],'fig');
    end
end

%% All DUSP1 distributions
if wantDUSP1DistPlots&&~isempty(modelChoice)
    % Create a figure
    % clear AllDuspDistsFig

    dexConc = {'0p3','1','10','100'};
    tTPLA = {'0','20','75','150','180'};

    for iDex = 1:9
        switch iDex
            case 1
                tDusp = [0    30    50    75    90   120   180];
                tSkip = 0;
            case 2
                tDusp = [0    30    50    75    90   120   180];
                tSkip = 0;
            case 3
                tDusp = [0    30    50    75    90   120   180];
                tSkip = 0;
            case 4
                tDusp = [0    10    20    30    40    50    60    75    90   120   150   180];
                tSkip = 0;
            case 5
                tDusp = [0    15    30    60];
                tSkip = 0;
            case 6
                tDusp = [0    10    20    35    50    80];
                tSkip = 2;
            case 7
                tDusp = [0    10    20    30    40    50    60    75    85    90   105   135];
                tSkip = 7;
            case 8
                tDusp = [0    10    20    30    40    50    60    75    90   120   150   160   180   210];
                tSkip = 10;
            case 9
                tDusp = [0    10    20    30    40    50    60    75    90   120   150   180   195   210   240];
                tSkip = 11;
        end
        colsMax = length(tDusp)-tSkip;

        AllDuspDistsFig{iDex} = figure('Visible',figVis);
        set(AllDuspDistsFig{iDex},'Position',[566         833        1359         522])

        % Define the number of rows and columns
        rows = 3;
        cols = max(12,colsMax);

        leftMargin = 0.08;
        bottomMargin = 0.05;
        topMargin = 0.02;
        width = 0.86*(1-leftMargin)/cols;
        colsep = 0.08*(1-leftMargin)/cols;
        heightB = 0.8*(1-bottomMargin-topMargin)/rows;
        rowsep = 0.18*(1-bottomMargin)/rows;
        fontsize = 13;

        rowLab = {'NucDist','CytDist','JointNucCyt'};

        rowLabels = {'Nucleus','Cytoplasm','Cytoplasm'};
        % Create axes in a loop
        for r = 1:3

            if iDex<=4
                oldFig = openfig(['savedFigures/RawFigs/',modelChoice,'/Dusp1_SSA_',dexConc{iDex},'_',rowLab{r},'.fig'],openFigVis);
            else
                oldFig = openfig(['savedFigures/RawFigs/',modelChoice,'/Dusp1_SSA_TPL_',tTPLA{iDex-4},'min_',rowLab{r},'.fig'],openFigVis);
            end

            for c = 1:colsMax

                % figure(AllDuspDistsFig{iDex})
                set(groot,'CurrentFigure',AllDuspDistsFig{iDex});

                % Calculate position for each axes
                axs{r,c} = axes('Position', [(c-1)*(width+colsep)+leftMargin, (rows-r)*(rowsep+heightB)+bottomMargin, width, heightB]);

                if iDex<=4
                    copyobj(oldFig.Children(length(oldFig.Children)-c+1).Children, axs{r,c});
                else
                    copyobj(oldFig.Children(length(oldFig.Children)-c+1-tSkip).Children, axs{r,c});
                end

                switch r
                    case 1
                        set(gca,'YLim',[0,0.05],'XLim',[-15,280]);
                    case 2
                        set(gca,'YLim',[0,0.03],'XLim',[-15,280]);
                    case 3
                        set(gca,'XLim',[-15,280],'YLim',[-15,280]);
                end

                if c~=1
                    set(gca,'yticklabels',[]);
                else
                    ylab = ylabel(rowLabels{r});
                    xlim = get(gca,'XLim');
                    ylab.Position(1) = -0.4*(xlim(2)-xlim(1));
                end
                %
                if r==1
                    title(['t = ',num2str(tDusp(c+tSkip)),' min']);
                end
                %
                % % % Minimize space between axes
                set(axs{r,c}, 'Box', 'on', 'YGrid', 'on', 'XGrid', 'on'); % Optional: add box around axes
            end
        end
    end
    % Adjust figure properties to minimize space
    set(gcf, 'Color', 'w'); % Set figure background color
    saveas(AllDuspDistsFig{1},['savedFigures/',modelChoice,'/Dusp1Dists_0p3'],'fig');
    saveas(AllDuspDistsFig{2},['savedFigures/',modelChoice,'/Dusp1Dists_1'],'fig');
    saveas(AllDuspDistsFig{3},['savedFigures/',modelChoice,'/Dusp1Dists_10'],'fig');
    saveas(AllDuspDistsFig{4},['savedFigures/',modelChoice,'/Dusp1Dists_100'],'fig');
    saveas(AllDuspDistsFig{5},['savedFigures/',modelChoice,'/Dusp1Dists_TPL_0'],'fig');
    saveas(AllDuspDistsFig{6},['savedFigures/',modelChoice,'/Dusp1Dists_TPL_20'],'fig');
    saveas(AllDuspDistsFig{7},['savedFigures/',modelChoice,'/Dusp1Dists_TPL_75'],'fig');
    saveas(AllDuspDistsFig{8},['savedFigures/',modelChoice,'/Dusp1Dists_TPL_150'],'fig');
    saveas(AllDuspDistsFig{9},['savedFigures/',modelChoice,'/Dusp1Dists_TPL_180'],'fig');


    out.figures.AllDuspDistsFig = AllDuspDistsFig;

    f = openfig(['savedFigures/',modelChoice,'/Dusp1Dists_100.fig']);
    axOld = findall(f,'Type','axes');
    axOldSort = axOld(end:-1:1);
    vv = [1,4,6,8,9,10,12];
    vvv=[vv,12+vv,24+vv];
    axOldSubset = axOldSort(vvv);
    newFig = figure;
    copyobj(axOldSubset,newFig);
    for i = 1:7
        set(newFig.Children(i),'Position',axOldSort(i).Position);
        set(newFig.Children(i+7),'Position',axOldSort(i+12).Position);
        set(newFig.Children(i+14),'Position',axOldSort(i+24).Position);
        newFig.Children(i).Color = [.91 .74 .85];%(1-C.CytoDUSP1_data)*0.15 + 0.85;
        newFig.Children(i+7).Color = [.83 .74 .89];%(1-C.NucDUSP1_data)*0.15 + 0.85;
        newFig.Children(i+14).Color = [.94 .81 .7];%(1-.5*C.NucDUSP1_data-.5*C.CytoDUSP1_data)*0.15 + 0.85;
        newFig.Children(i).Children(2).EdgeColor = C.NucDUSP1_data;
        newFig.Children(i+7).Children(2).EdgeColor = C.CytoDUSP1_data;

    end
    for i = 1:length(newFig.Children)
        newFig.Children(i).Position(1,3) = newFig.Children(i).Position(1,3)*12/7;
        newFig.Children(i).Position(1,1) = newFig.Children(1).Position(1,1) + (newFig.Children(i).Position(1,1)-newFig.Children(1).Position(1,1))*12/7;
        newFig.Children(i).FontSize = smallFont;
        newFig.Children(i).FontWeight = 'bold';
        newFig.Children(i).XLabel.Visible = 'off';
        newFig.Children(i).YLabel.Visible = 'off';
    end

    newFig.Position = [1018         575         706         472];
    if strcmp(modelChoice,finalFullModel)
        saveas(newFig,'savedFigures/FinalFigures/MainFigures/Fig3C_DUSP1Distributions.fig','fig');
    end
    saveas(f,['savedFigures/FinalFigures/SupplementalFigures/Fit_DUSP1Distributions',modelChoice,'.fig'],'fig');


    %% SI DUSP1 Prediction Distribution Plots
    bigFig = figure; clf; figure(bigFig)
    set(bigFig,'Position',[19         558        1439         741]);
    k = 0;
    for i = 1:4
        copyobj(AllDuspDistsFig{i}.Children,bigFig)
        for j = 1:length(AllDuspDistsFig{i}.Children)
            bigFig.Children(k+j).Position(4) = bigFig.Children(k+j).Position(4)*0.23;
            bigFig.Children(k+j).Position(2) = 0.03+0.23*(4-i)+bigFig.Children(k+j).Position(2)/4;
        end
    end
    for i = 1:length(bigFig.Children)
        bigFig.Children(i).FontSize = smallFont;
        bigFig.Children(i).FontWeight = 'bold';
        if i<length(bigFig.Children)-12
            bigFig.Children(i).Title.Visible = 'off';
        end
        if i>12
            bigFig.Children(i).XTickLabel = [];
        end
    end
    for i = [26,30,32,34,35]
        bigFig.Children(i).Title.Visible = 'on';
    end
    vShift = [1 3 4 5 7 9 12];
    for j = 1:9
        for i = 1:7
            bigFig.Children(36+(j-1)*7+i).Position(1) = bigFig.Children(vShift(i)).Position(1);
        end
    end
    for i =[1:12,[44:50,65:71,86:92]-7]
        bigFig.Children(i).YLim = [0,200];
        bigFig.Children(i).XLim = [0,150];
    end
    for i =[13:24,[44:50,65:71,86:92]]
        bigFig.Children(i).YLim = [0,0.02];
        bigFig.Children(i).XLim = [0,150];
    end
    for i =[25:36,[44:50,65:71,86:92]+7]
        bigFig.Children(i).YLim = [0,0.04];
        bigFig.Children(i).XLim = [0,150];
    end
    for i =[12,36+7,36+7+21,36+7+42]
        bigFig.Children(i).YLabel.String = {'Joint','Cyto (y) vs.','Nuc (x)'};
        bigFig.Children(i).YLabel.Position(2) = 20;
        bigFig.Children(i).YLabel.Position(1) = -100;
        bigFig.Children(i).YLabel.Rotation = 0;
    end
    for i =[24,36+14,36+14+21,36+14+42]
        bigFig.Children(i).YLabel.String = 'Cytoplasm';
        bigFig.Children(i).YLabel.Position(2) = 0.007;
        bigFig.Children(i).YLabel.Position(1) = -100;
        bigFig.Children(i).YLabel.Rotation = 0;
    end
    for i =[36,36+21,36+21+21,36+21+42]
        bigFig.Children(i).YLabel.String = 'Nucleus';
        bigFig.Children(i).YLabel.Position(2) = 0.016;
        bigFig.Children(i).YLabel.Position(1) = -100;
        bigFig.Children(i).YLabel.Rotation = 0;
    end
    bigFig.Children(7).XLabel.String = '       - - - - - - - - - - - - - - - - - - - - - - - - - - - - Number of Molecules - - - - - - - - - - - - - - - - - - - - - - - - - - - - ';
    
    saveas(f,['savedFigures/FinalFigures/SupplementalFigures/All_DUSP1Distributions',modelChoice,'.fig'],'fig');
        
    %%
end

%% Make MH plots for GR Model
if wantMHGR
    close all
    [combined_GRModel,log10PriorMean,log10PriorStd] = dusp1ModelLibrary_Final('combinedGRModel');
    load('savedParameters/MHResults_GR_Model.mat')
    names = {'$k_{cn0}$'
        '$k_{cn1}$'
        '$g_{Dex}$'
        '$k_{nc}$'
        '$k_{GR}$'
        '$\gamma_{GRcyt}$'
        '$\gamma_{GRnuc}$'
        '$M_{Dex}$'};
    descriptions = {'GR transport Cyt$\\rightarrow$Nuc, basal (min$^{-1}$)',...
'GR transport Cyt$\\rightarrow$Nuc, Dex-dependent (min$^{-1}$)',...
'Dex degradation rate, effectively zero (min$^{-1}$)',...
'GR transport Nuc$\\rightarrow$Cyt (min$^{-1}$)',...
'GR synthesis rate, Cyt (min$^{-1}$)',...
'GR degradation, Cyt (min$^{-1}$)',...
'GR degradation, Nuc (min$^{-1}$)' ,...
'Michaelis constant for Dex-induced import (nM)'};
    combined_GRModel.SSITModels{1}.plotMHResults(mhResultsGR,ESS=false,names=names,...
        latexFileName='ParameterTables/MHResultTableGR.txt',descriptions=descriptions,...
        showMarginalPosteriors=true,...
        priorMean=log10PriorMean(combined_GRModel.SSITModels{1}.fittingOptions.modelVarsToFit),...
        priorSig=log10PriorStd(combined_GRModel.SSITModels{1}.fittingOptions.modelVarsToFit));

    f=figure(4);saveas(f,'savedFigures/FinalFigures/SupplementalFigures/MarginalPosteriorDistributions_GR.fig','fig')
    f=figure(3);saveas(f,'savedFigures/FinalFigures/SupplementalFigures/JointPosteriorDistributions_GR.fig','fig')
end

%% Make MH Plots for Full Model - Semi-Mechanistic
if wantMHCyt
    close all
    [dusp1Model,log10PriorMean,log10PriorStd] = dusp1ModelLibrary_Final('ModelDUSP1_100nM');
    dusp1Model.fittingOptions.modelVarsToFit = [1,2,3,4,14,15,16,19,20,21,22,24,25,26];
    % try 
        load('savedParameters/MHResults_FinalModel_allButGR.mat')
    % catch
    %     load TMPMHFinalModel_allButGR
    %     MHResultsAll.mhSamples = smpl(value~=0,:);
    %     MHResultsAll.mhValue = value(value~=0);
    % end
    names = {'$k_{off}$'
        '$k_{on,0}$'
        '$kr_{on,0}$'
        '$k_{nc}$'
        '$k_{degCyt,0}$' %5
        '$k_{degCyt,1}$'
        '$k_{DTTP}$'
        '$\tau_{elong}$'
        '$TTP_1$'
        '$t_{TTP}$' %10
        '$\eta$'
        '$k_{on,1}$'
        '$m_{koff}$'
        '$k_{r,on,1}$' };
    descriptions = {'Basal gene inactivation rate (min$^{-1}$)' ;
    'Basal gene activation rate (min$^{-1}$)'              ;
    'nucGR-dependent gene activation rate (min$^{-1}$)'    ;
    'nucGR modulation of inactivation (molecules$^{-1}$)'  ;
    'Basal DUSP1 transcription rate (min$^{-1}$)'          ;
    'nucGR-dependent DUSP1 transcription rate (min$^{-1}$)';
    'DUSP1 elongation time (min)'                         ;
    'DUSP1 nuclear export rate (min$^{-1}$)'               ;
    'Basal DUSP1 mRNA degradation rate (min$^{-1}$)'       ;
    'Max TTP-dependent degradation rate (min$^{-1}$)'      ;
    'Max fold-change in TTP-mediated degradation'         ;
    'Half-max time for TTP activity (min)'                ;
    'Hill coefficient, TTP time-dependence'               ;
    'TTP-mRNA binding saturation (molecules)'              };

    parameterReorder = [1,2,12,13,3,14,8,4,5,6,9,10,11,7];

    dusp1Model.plotMHResults(MHResultsAll,ESS=false,names=names,...
        latexFileName='ParameterTables/MHResultTableSemiMech.txt',...
        parameterReorder = parameterReorder, descriptions=descriptions,...
        showMarginalPosteriors=true,...
        priorMean=log10PriorMean(dusp1Model.fittingOptions.modelVarsToFit),...
        priorSig=log10PriorStd(dusp1Model.fittingOptions.modelVarsToFit));
    f=figure(5);saveas(f,'savedFigures/FinalFigures/SupplementalFigures/MarginalPosteriorDistributions_SemiMech.fig','fig')
    f=figure(3);saveas(f,'savedFigures/FinalFigures/SupplementalFigures/JointPosteriorDistributions_SemiMech_1to49.fig','fig')
    f=figure(4);saveas(f,'savedFigures/FinalFigures/SupplementalFigures/JointPosteriorDistributions_SemiMech_50to91.fig','fig')
end

%%
if wantMHMech
    close all
    mechModel = dusp1ModelLibrary_Final('GR_DUSP1_TTP_Model_2');
    mechModel.fittingOptions.modelVarsToFit = [1,2,3,4,14,15,19,24,25,26,27:36];
    % try
        load('savedParameters/MHResults_Pars_GR_DUSP1_TTP_Model_All.mat')
    % catch
    %     load TMPMHPars_GR_DUSP1_TTP_Model_All
    %     MHResultsAll.mhSamples = smpl(value~=0,:);
    %     MHResultsAll.mhValue = value(value~=0);
    % end
    names = {'$k_{off}$'
        '$k_{on,0}$'
        '$kr_{on,0}$'
        '$k_{nc}$'
        '$k_{degCyt,0}$' %5
        '$k_{degCyt,1}$'
        '$\tau_{elong}$'
        '$TTP_1$'
        '$t_{TTP}$'
        '$\eta$'
        '$k_{on,1}$'
        '$m_{koff}$'
        '$k_{r,on,1}$' 
        'kon_{TTP,1}'   
        'koff_{TTP}'    
        'kr_{on,TTP}'   
        'k_{TTP,nc}'
        'k_{translate}'  
        '\gamma_{TTP-Prot}'     
        'k_{bind}'
        'eta_{TTP}'     
        'k_{unbind}'      
        '\gamma_{TTP-RNA}'};

    descriptions = {'Basal gene inactivation rate (min$^{-1}$)','$k_{off}$'        ;
        'Basal gene activation rate (min$^{-1}$)'              ,'$k_{on,0}$';
        'nucGR-dependent gene activation rate (min$^{-1}$)'    ,'$k_{on,1}$';
        'nucGR modulation of inactivation (molecules$^{-1}$)'  ,'$m_{koff}$';
        'Basal DUSP1 transcription rate (min$^{-1}$)'          ,'$kr_{on,0}$';
        'nucGR-dependent DUSP1 transcription rate (min$^{-1}$)','$kr_{on,1}$';
        'DUSP1 elongation time (min)'                          ,'$\tau_{elong}$';
        'DUSP1 nuclear export rate (min$^{-1}$)'               ,'$k_{nc}$';
        'Basal DUSP1 mRNA degradation rate (min$^{-1}$)'       ,'$\gamma_{Cyt,0}$';
        'Max TTP-dependent degradation rate (min$^{-1}$)'      ,'$\gamma_{Cyt,1}$';
        % 'Max fold-change in TTP-mediated degradation'          ,'$TTP1$';
        % 'Half-max time for TTP activity (min)'                 ,'$t_{TTP}$';
        % 'Hill coefficient, TTP time-dependence'                ,'$\eta$';
        % 'TTP-mRNA binding saturation (molecules)'              ,'$kD_{TTP}$';
        'nucGR-dependent TTP gene activation rate (min$^{-1}$)','$kon_{TTP,1}$';
        'TTP gene inactivation rate (min$^{-1}$)'              ,'$koff_{TTP}$';
        'TTP transcription rate (min$^{-1}$)'                  ,'$kr_{on,TTP}$';
        'TTP nuclear export rate (min$^{-1}$)'                 ,'$k_{TTP,nc}$';
        'Basal TTP mRNA degradation rate (min$^{-1}$)'         ,'$\gamma_{TTP-RNA}$';
        'TTP translation rate (min$^{-1}$)'                    ,'$k_{translate}$';
        'TTP protein degradation rate (min$^{-1}$)'            ,'$\gamma_{TTP-Prot}$';
        'TTP–mRNA binding rate (min$^{-1}$ mol$^{-1}$)'        ,'$k_{bind}$';
        'TTP–mRNA unbinding rate (min$^{-1}$)'                 ,'$k_{unbind}$';
        'Hill coefficient, TTP–mRNA binding'                   ,'$\eta_{TTP}$'     };

        parameterReorder = [1,2,8,9,3,10,7,4,5,6,11,12,13,14,20,15,16,17,19,18];

    mechModel.plotMHResults(MHResultsAll,ESS=false,names=names,...
        descriptions = descriptions, parameterReorder=parameterReorder, ...
        latexFileName='ParameterTables/MHResultTableMechanistic.txt',...
        showMarginalPosteriors=true,...
        priorMean=log10PriorMean(mechModel.fittingOptions.modelVarsToFit),...
        priorSig=log10PriorStd(mechModel.fittingOptions.modelVarsToFit));
    f=figure(7);saveas(f,'savedFigures/FinalFigures/SupplementalFigures/MarginalPosteriorDistributions_Mech.fig','fig')
    f=figure(3);saveas(f,'savedFigures/FinalFigures/SupplementalFigures/JointPosteriorDistributions_Mech_1to49.fig','fig')
    f=figure(4);saveas(f,'savedFigures/FinalFigures/SupplementalFigures/JointPosteriorDistributions_Mech_50to98.fig','fig')
    f=figure(5);saveas(f,'savedFigures/FinalFigures/SupplementalFigures/JointPosteriorDistributions_Mech_99to147.fig','fig')
    f=figure(6);saveas(f,'savedFigures/FinalFigures/SupplementalFigures/JointPosteriorDistributions_Mech_148to190.fig','fig')
end

%% LOCAL FUNCTIONS

end

function selectedPlots = normalizePlotSelection(plots)
allPlots = {'GRdistributions','GRmeans','TotalGR','GRModelComparison', ...
    'GRmeansAndVars','triptolide','allGRdistributions','DUSP1DistPlots',...
    'MHGR','MHCyt','MHMech'};

if ischar(plots)
    plots = string({plots});
elseif isstring(plots)
    plots = reshape(plots, 1, []);
elseif iscell(plots)
    plots = string(plots);
else
    error('opts.plots must be text or a collection of text values.');
end

plots = lower(strrep(strtrim(plots), ' ', ''));
selectedPlots = {};

for i = 1:numel(plots)
    key = char(plots(i));
    switch key
        case 'all'
            selectedPlots = allPlots;
            return
        case {'grdistributions','grdist','grdists'}
            selectedPlots{end+1} = 'GRdistributions';
        case {'grmeans','grmean'}
            selectedPlots{end+1} = 'GRmeans';
        case {'grmeansandvars','grmeansandvariances','grmeansandvar','grmeansandvariance'}
            selectedPlots{end+1} = 'GRmeansAndVars';
        case {'triptolidemeasurements','tplmeasurements','tplmsmts','triptolide'}
            selectedPlots{end+1} = 'triptolide';
        case {'allgrdistributions','allgrdists','allgr'}
            selectedPlots{end+1} = 'allGR';
        case {'DUSP1DistPlotsdistributions','dusp1distributions','DUSP1DistPlotsdists','dusp1dists','dusp1distplots'}
            selectedPlots{end+1} = 'DUSP1DistPlots';
        case {'totalgrmeans','totalgr','grtotal','totalgrmean'}
            selectedPlots{end+1} = 'TotalGR';
        case {'grmodelcomparison','grmodelcompare','comparegrmodels','grcomparison'}
            selectedPlots{end+1} = 'GRmodelComparison';
        case {'mhgr'}
            selectedPlots{end+1} = 'MHGR';
        case {'mhcyt'}
            selectedPlots{end+1} = 'MHCyt';
        case {'mhmech'}
            selectedPlots{end+1} = 'MHMech';

        otherwise
            error('Unknown plot selection: %s', key);
    end
end

selectedPlots = unique(selectedPlots, 'stable');
end

function modelChoice = normalizeFinalModelChoice(modelChoice, savedFullParSet)
if isnumeric(modelChoice)
    switch modelChoice
        case {2}
            modelChoice = savedFullParSet{1};
        case {10}
            modelChoice = savedFullParSet{2};
        case {3}
            modelChoice = savedFullParSet{3};
        case {4}
            modelChoice = savedFullParSet{4};
        case {8}
            modelChoice = savedFullParSet{5};
        otherwise
            modelChoice=[];%error('Unknown modelChoice index: %g', modelChoice);
    end
else
    error('Function disabled')
    % modelChoice = char(string(modelChoice));
    % aliases = containers.Map( ...
    %     {'kon','kon_only','simple','full','allthree','mechanistic','ttp'}, ...
    %     {savedFullParSet{1}, savedFullParSet{1}, savedFullParSet{1}, ...
    %      savedFullParSet{2}, savedFullParSet{2}, savedFullParSet{3}, savedFullParSet{3}});
    % key = lower(strrep(modelChoice,' ','_'));
    % if isKey(aliases,key)
    %     modelChoice = aliases(key);
    % end
    % if ~ismember(modelChoice, savedFullParSet)
    %     error('Unknown modelChoice: %s', modelChoice);
    % end
end
end

function ax = formatLines(ax,col)
% Standardize line widths and errorbar markers.
% If 'col' is provided (char color or [r g b]), color ONLY data errorbars
% and tint any patches; leave model mean 'line' objects to be set elsewhere.
arguments
    ax
    col = []
end

% Decide if we will color with numeric or char
useNumeric = isnumeric(col) && ~isempty(col);

% Iterate children (top to bottom)
for i = 1:length(ax.Children)
    t = get(ax.Children(i), 'Type');
    switch t
        case 'errorbar'
            ax.Children(i).LineWidth = 2;
            if isprop(ax.Children(i),'MarkerSize')
                ax.Children(i).MarkerSize = 8;
            end
            if isprop(ax.Children(i),'Marker')
                mk = ax.Children(i).Marker;
            else
                mk = '';
            end
            % Only set marker if not already none/empty
            if ~(~isempty(mk) && strcmp(mk,'none'))
                ax.Children(i).Marker = 'o';
            end

            if ~isempty(col)
                if useNumeric
                    if ~strcmp(mk,'none')
                        ax.Children(i).Color = col;
                        if isprop(ax.Children(i),'MarkerFaceColor'), ax.Children(i).MarkerFaceColor = col; end
                        if isprop(ax.Children(i),'MarkerEdgeColor'), ax.Children(i).MarkerEdgeColor = col; end
                    end
                else
                    if ~strcmp(mk,'none')
                        ax.Children(i).Color = col;
                        if isprop(ax.Children(i),'MarkerFaceColor'), ax.Children(i).MarkerFaceColor = col; end
                        if isprop(ax.Children(i),'MarkerEdgeColor'), ax.Children(i).MarkerEdgeColor = col; end
                    end
                end
            end

        case 'line'
            ax.Children(i).LineWidth = 2; % do NOT recolor here

        case 'patch'
            if ~isempty(col)
                if useNumeric
                    blend = 0.85; % local default; separate from C.Shade_blend
                    face  = blend*ones(1,3) + (1-blend)*col;
                    if isprop(ax.Children(i),'FaceColor'), ax.Children(i).FaceColor = face; end
                else
                    switch col
                        case 'r'
                            ax.Children(i).FaceColor = [1,0.9,0.9];
                        case 'b'
                            ax.Children(i).FaceColor = [0.9,0.9,1];
                        otherwise
                            % leave as-is
                    end
                end
            end
    end
end

end


function recolorAxes(ax, args)
% recolorAxes(ax, struct(...))
% Fields:
%   dataColor      : [r g b] for data errorbars/marker lines (with markers)
%   modelLine      : [r g b] for model mean line(s) (Marker == 'none')
%   shadeFromModel : true/false (if true, tint patches from modelLine color)
%   shadeBlend     : scalar in [0,1] for how light the patch should be

if ~ishandle(ax) || ~strcmp(get(ax,'Type'),'axes')
    return
end

ch = ax.Children;

% DATA ERRORBARS + DATA LINES WITH MARKERS
if isfield(args,'dataColor') && ~isempty(args.dataColor)
    eb = findobj(ch,'Type','errorbar');
    for k = 1:numel(eb)
        mk = '';
        if isprop(eb(k),'Marker'), mk = eb(k).Marker; end
        if ~strcmp(mk,'none') % treat only marker-bearing errorbars as data
            eb(k).Color = args.dataColor;
            if isprop(eb(k),'MarkerFaceColor'), eb(k).MarkerFaceColor = args.dataColor; end
            if isprop(eb(k),'MarkerEdgeColor'), eb(k).MarkerEdgeColor = args.dataColor; end
        end
    end
    datalines = findobj(ch,'Type','line','-not','Marker','none');
    for k = 1:numel(datalines)
        datalines(k).Color = args.dataColor;
    end
end

% MODEL MEAN LINES (no markers)
if isfield(args,'modelLine') && ~isempty(args.modelLine)
    modellines = findobj(ch,'Type','line','-and','Marker','none');
    for k = 1:numel(modellines)
        modellines(k).Color = args.modelLine;
    end
end

% MODEL SHADING (patch)
if isfield(args,'shadeFromModel') && args.shadeFromModel
    patches = findobj(ch,'Type','patch');
    if ~isempty(patches)
        base = [0 0 0];
        if isfield(args,'modelLine') && ~isempty(args.modelLine), base = args.modelLine; end
        blend = 0.85;
        if isfield(args,'shadeBlend') && ~isempty(args.shadeBlend), blend = args.shadeBlend; end
        face = blend*ones(1,3) + (1-blend)*base;
        for k = 1:numel(patches)
            if isprop(patches(k),'FaceColor'), patches(k).FaceColor = face; end
        end
    end
end
end


function listChildren(ax)
% listChildren(ax)
% Utility to print index/type/visibility/marker/name/colors for each child.
ch = ax.Children;
N = numel(ch);
fprintf('Children for axes (top to bottom = 1..%d):\n', N);
for k = 1:N
    t = get(ch(k),'Type');
    vis = get(ch(k),'Visible','on');
    dn  = '';
    mk  = '';
    col = '';
    fcol= '';
    if isprop(ch(k),'DisplayName'), dn = ch(k).DisplayName; end
    if isprop(ch(k),'Marker'),      mk = ch(k).Marker;      end
    if isprop(ch(k),'Color'),       col = mat2str(ch(k).Color,3); end
    if isprop(ch(k),'FaceColor'),   fcol= mat2str(ch(k).FaceColor,3); end
    fprintf('%2d) %-9s  Vis:%-3s  Marker:%-4s  Name:%s  Color:%s  Face:%s\n',...
        k, t, vis, mk, dn, col, fcol);
end
end

function ax = plotTotalGRFromMeanFigures(ax, fNuc, fCyt, ratio_N2C, C)
% plotTotalGRFromMeanFigures
% Combines nuclear and cytoplasmic GR mean/variance plots into Total GR.
%
% The GR mean figures show nuclear and cytoplasmic GR as concentrations.
% To combine them into total mass-like GR, nuclear GR is multiplied by the
% nucleus-to-cytoplasm ratio before adding to cytoplasmic GR.

if isempty(fNuc) || isempty(fCyt) || ~ishandle(fNuc) || ~ishandle(fCyt)
    error('Valid nuclear and cytoplasmic GR figure handles are required.');
end

if ~ishandle(ax) || ~strcmp(get(ax, 'Type'), 'axes')
    error('A valid target axes handle is required.');
end

axNuc = findMainAxes(fNuc);
axCyt = findMainAxes(fCyt);

% Find data objects robustly.
hNucErr = findobj(axNuc, 'Type', 'ErrorBar', 'Visible', 'on');
hCytErr = findobj(axCyt, 'Type', 'ErrorBar', 'Visible', 'on');

if isempty(hNucErr) || isempty(hCytErr)
    error('Could not find nuclear/cytoplasmic errorbar objects for Total GR plot.');
end

% Prefer marker-bearing errorbars, since those are the experimental data.
hNucErr = pickDataErrorbar(hNucErr);
hCytErr = pickDataErrorbar(hCytErr);

xNuc    = hNucErr.XData;
yNuc    = hNucErr.YData;
ynegNuc = hNucErr.YNegativeDelta;
yposNuc = hNucErr.YPositiveDelta;

xCyt    = hCytErr.XData;
yCyt    = hCytErr.YData;
ynegCyt = hCytErr.YNegativeDelta;
yposCyt = hCytErr.YPositiveDelta;

% Find model uncertainty bands robustly.
pNuc = findobj(axNuc, 'Type', 'Patch', 'Visible', 'on');
pCyt = findobj(axCyt, 'Type', 'Patch', 'Visible', 'on');

if isempty(pNuc) || isempty(pCyt)
    error('Could not find nuclear/cytoplasmic patch objects for Total GR plot.');
end

pNuc = pNuc(1);
pCyt = pCyt(1);

xPolyNuc = pNuc.XData(:)';
yPolyNuc = pNuc.YData(:)';
xPolyCyt = pCyt.XData(:)';
yPolyCyt = pCyt.YData(:)';

nNuc = numel(xPolyNuc)/2;
nCyt = numel(xPolyCyt)/2;

if mod(numel(xPolyNuc),2) ~= 0 || mod(numel(xPolyCyt),2) ~= 0
    error('Patch XData lengths must be even to infer lower/upper uncertainty bands.');
end

xModelNuc = xPolyNuc(1:nNuc);
xModelCyt = xPolyCyt(1:nCyt);

lowerNuc = yPolyNuc(1:nNuc);
upperNuc = yPolyNuc(end:-1:nNuc+1);

lowerCyt = yPolyCyt(1:nCyt);
upperCyt = yPolyCyt(end:-1:nCyt+1);

% Convert nuclear GR to cytoplasm-scaled total amount before combining.
muNuc  = (lowerNuc + upperNuc)/2 * ratio_N2C;
stdNuc = (upperNuc - lowerNuc)/2 * ratio_N2C;

muCyt  = (lowerCyt + upperCyt)/2;
stdCyt = (upperCyt - lowerCyt)/2;

if numel(muNuc) ~= numel(muCyt)
    error('Nuclear and cytoplasmic model bands have different lengths.');
end

muTot  = muNuc + muCyt;
stdTot = sqrt(stdNuc.^2 + stdCyt.^2);

xModel = xModelCyt;

% Combine experimental data.
if numel(yNuc) ~= numel(yCyt)
    error('Nuclear and cytoplasmic data vectors have different lengths.');
end

xData   = xCyt;
yData   = yCyt + ratio_N2C*yNuc;
ynegTot = sqrt(ynegCyt.^2 + (ratio_N2C*ynegNuc).^2);
yposTot = sqrt(yposCyt.^2 + (ratio_N2C*yposNuc).^2);

% axes(ax);
cla(ax);
hold(ax, 'on');

fill(ax, ...
    [xModel, xModel(end:-1:1)], ...
    [muTot-stdTot, muTot(end:-1:1)+stdTot(end:-1:1)], ...
    C.Model_line, ...
    'FaceAlpha', 0.25, ...
    'EdgeColor', 'none');

errorbar(ax, xData, yData, ynegTot, yposTot, 'o', ...
    'Color', C.TotalGR_data, ...
    'MarkerFaceColor', C.TotalGR_data, ...
    'MarkerEdgeColor', C.TotalGR_data, ...
    'LineWidth', 1.2);

plot(ax, xModel, muTot, ...
    'Color', C.Model_line, ...
    'LineWidth', 2.5);

ax = formatLines(ax);
hold(ax, 'off');

end

function ax = findMainAxes(figHandle)
% Find the main plotting axes while avoiding legends/colorbars when possible.

axs = findobj(figHandle, 'Type', 'Axes');

if isempty(axs)
    error('No axes found in the supplied figure.');
end

% Prefer axes that actually contain plotted graphics.
for k = 1:numel(axs)
    ch = axs(k).Children;
    if any(arrayfun(@(h) any(strcmp(get(h,'Type'), ...
            {'line','errorbar','patch'})), ch))
        ax = axs(k);
        return
    end
end

% Fallback.
ax = axs(1);
end

function h = pickDataErrorbar(hErr)
% Prefer marker-bearing errorbar objects as experimental data.

h = hErr(1);

for k = 1:numel(hErr)
    mk = '';
    if isprop(hErr(k), 'Marker')
        mk = hErr(k).Marker;
    end

    if ~isempty(mk) && ~strcmp(mk, 'none')
        h = hErr(k);
        return
    end
end
end

function [dataOut, modelOut] = extractGRMeanVarFromFig(figHandle, compartment, ratio_N2C)
% extractGRMeanVarFromFig
%
% Extracts experimental data errorbars and model mean line from one of the
% saved GR mean/variance figures, e.g.,
%   combinedGRModel_Dex100_3.fig
%
% compartment:
%   'cyt' or 'nuc'

ax = findMainAxes(figHandle);

ch = ax.Children;

switch lower(compartment)
    case {'nuc','nuclear'}
        Plts = [1,3,4,2,5,6];
        ratio = ratio_N2C;
    case {'cyt','cyto','cytoplasmic'}
        Plts = [2,5,6,1,3,4];
        ratio = 1;
    otherwise
        error('Unknown compartment "%s". Use "nuc" or "cyt".', compartment);
end

% Reorder exactly as the existing makeFinalFigures GR mean helper does
ch = ch(Plts);
ch = ch(1:min(3,numel(ch)));

% In the existing code, after this reordering:
%   ch(1) is the data errorbar
%   ch(2) and/or ch(3) are model-related objects
% depending on MATLAB object ordering. We extract robustly below.

% Experimental data: marker-bearing errorbar
errObjs = findobj(ch, 'Type', 'ErrorBar');
if isempty(errObjs)
    error('Could not find experimental errorbar object.');
end
hData = pickDataErrorbar(errObjs);

dataOut.x = hData.XData(:)';
dataOut.y = hData.YData(:)' ./ ratio;

if isprop(hData, 'YNegativeDelta') && ~isempty(hData.YNegativeDelta)
    dataOut.yneg = hData.YNegativeDelta(:)' ./ ratio;
else
    dataOut.yneg = zeros(size(dataOut.y));
end

if isprop(hData, 'YPositiveDelta') && ~isempty(hData.YPositiveDelta)
    dataOut.ypos = hData.YPositiveDelta(:)' ./ ratio;
else
    dataOut.ypos = zeros(size(dataOut.y));
end

% Model mean: prefer a line with Marker == 'none'. If multiple exist,
% choose the one with the most XData points.
lineObjs = findobj(ch, 'Type', 'Line');

if isempty(lineObjs)
    error('Could not find model mean line.');
end

bestIdx = 1;
bestN = -inf;

for k = 1:numel(lineObjs)
    mk = '';
    if isprop(lineObjs(k), 'Marker')
        mk = lineObjs(k).Marker;
    end

    if strcmp(mk, 'none') || isempty(mk)
        nPts = numel(lineObjs(k).XData);
        if nPts > bestN
            bestN = nPts;
            bestIdx = k;
        end
    end
end

hModel = lineObjs(bestIdx);

modelOut.x = hModel.XData(:)';
modelOut.y = hModel.YData(:)' ./ ratio;

% Sort model curve by time:
[modelOut.x, sortIdx] = sort(modelOut.x);
modelOut.y = modelOut.y(sortIdx);

% Sort data by time:
[dataOut.x, sortIdx] = sort(dataOut.x);
dataOut.y = dataOut.y(sortIdx);
dataOut.yneg = dataOut.yneg(sortIdx);
dataOut.ypos = dataOut.ypos(sortIdx);

end
