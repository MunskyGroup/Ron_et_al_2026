%% DUSP1 Triptolide Cyto Degradation - data pulled from saved .fig
% clear
% figfileDusp1Tpl = 'savedFigures/FinalModel_allButGR_NoPrior/TplMsmtsVsTime.fig';
% figfileDusp1Tpl = 'savedFigures/Pars_GR_DUSP1_TTP_Model_All_2/TplMsmtsVsTime.fig';

%% Parameters
ttpl_list = [0, 20, 75, 150, 180];

%% Open figure invisibly and locate the bottom (DUSP1-Cyt) row of axes
fig_in  = openfig(figfileDusp1Tpl, 'invisible');
cleanup = onCleanup(@() close(fig_in));
allAxes = findall(fig_in, 'Type', 'axes');
cytoAxes  = allAxes([5:-1:1]);
nucAxes  = allAxes([10:-1:6]);

%% Process Model at Different Times.
compartments = {'Nuc','Cyt'};
for kCompartment = 1:2
    switch kCompartment
        case 2
            Axes = cytoAxes;
        case 1
            Axes = nucAxes;
    end

    for iTPL = 1:5
        % Load Axis of plot made for this TPL time, and extract the mean and
        % patch limits found for the model.
        ax   = Axes(iTPL);
        pa = findobj(ax, 'Type', 'Patch');
        mn = findobj(ax, 'Type', 'Line');

        mod(iTPL).t_abs_Mod = mn(2).XData;
        mod(iTPL).y_abs_Mod = mn(2).YData;
        mod(iTPL).t_abs_Mod_pa = pa.XData';
        mod(iTPL).patch_abs_Mod = pa.YData';

        % Convert patch limits to deltas from the mean/median values.
        mod(iTPL).eps_abs_Mod = [-mod(iTPL).patch_abs_Mod(1:end/2)+mod(iTPL).y_abs_Mod,...
            mod(iTPL).patch_abs_Mod(end/2+1:end)-mod(iTPL).y_abs_Mod(end:-1:1)];

        % Extract the means and error bars from the Data
        eb = findobj(ax, 'Type', 'errorbar');
        dat(iTPL).y_abs_Dat = eb.YData;
        dat(iTPL).eps_abs_Dat = [eb.YNegativeDelta;eb.YPositiveDelta];

    end

    % Recenter time data to zero at time of the TPL application, and remove
    % earlier data, which is repeated from fitting.
    for iTPL = 1:length(mod)
        mod(iTPL).t_abs_Mod = mod(iTPL).t_abs_Mod - ttpl_list(iTPL);
        J = mod(iTPL).t_abs_Mod>=0;
        mod(iTPL).t_abs_Mod = mod(iTPL).t_abs_Mod(J);
        mod(iTPL).y_abs_Mod = mod(iTPL).y_abs_Mod(J);
        dat(iTPL).y_abs_Dat = dat(iTPL).y_abs_Dat(J);
        dat(iTPL).eps_abs_Dat = dat(iTPL).eps_abs_Dat(:,J);

        mod(iTPL).t_abs_Mod_pa = mod(iTPL).t_abs_Mod_pa - ttpl_list(iTPL);
        J = mod(iTPL).t_abs_Mod_pa>=0;
        mod(iTPL).t_abs_Mod_pa = mod(iTPL).t_abs_Mod_pa(J);
        mod(iTPL).patch_abs_Mod = mod(iTPL).patch_abs_Mod(J);
        mod(iTPL).eps_abs_Mod = mod(iTPL).eps_abs_Mod(J);
    end


    % Estimate Time at which Free Decay Begins (after delay from TPL transport,
    % Elongtion, and Nuclear export)
    kTPL = 0.098321;
    tau_elong = 1.76192;
    ktransport = 0.088409; %Temporary values -- need to use actual values from model
    delay = log(2)/kTPL + tau_elong + log(2)/ktransport;

    % Normalize Y-Axis by Dividing by Model Mean at Delay Time.
    for iTPL = 1:5
        % Extrapolate using times after delay to find mean of model at delay
        % time
        X = mod(iTPL).t_abs_Mod;
        J = X>=delay;
        X = X(J);
        Y = log(mod(iTPL).y_abs_Mod);
        Y = Y(J);
        Ydelay = exp(interp1(X,Y,delay,"linear","extrap"));

        % Apply scaling
        mod(iTPL).y_norm_Mod = mod(iTPL).y_abs_Mod / Ydelay;
        mod(iTPL).eps_norm_Mod = mod(iTPL).eps_abs_Mod / Ydelay;

        % Reform Normalized Patch Coordinates
        mod(iTPL).patch_norm_Mod = [mod(iTPL).y_norm_Mod,mod(iTPL).y_norm_Mod(end:-1:1)] +...
            [-mod(iTPL).eps_norm_Mod(1:end/2),mod(iTPL).eps_norm_Mod(end/2+1:end)];

        % Normalize the Data by Model Prediction at Delay Time
        dat(iTPL).y_norm_Dat = dat(iTPL).y_abs_Dat / Ydelay;
        dat(iTPL).eps_norm_Dat = dat(iTPL).eps_abs_Dat / Ydelay;
    end


    % Make Plot for Model in Each Experiment
    f = figure(100+10*kCompartment);clf
    set(f,'Name',[compartments{kCompartment},' DUSP1 vs Time after Triptolide']);
    for iTPL = 1:length(mod)
        % Original Model
        subplot(2,length(mod),iTPL)
        patch(mod(iTPL).t_abs_Mod_pa,mod(iTPL).patch_abs_Mod,'g'); hold on
        plot(mod(iTPL).t_abs_Mod,mod(iTPL).y_abs_Mod,'k','linewidth',2)

        % Add Data Error Bars Unnormalized
        errorbar(mod(iTPL).t_abs_Mod,dat(iTPL).y_abs_Dat,...
            dat(iTPL).eps_abs_Dat(1,:),dat(iTPL).eps_abs_Dat(2,:),'linestyle','--','LineWidth',2);

        set(gca,'YScale','lin','ylim',[0,200],'FontSize',16,'xtick',[0,20,40])

        % Normalized Model
        subplot(2,length(mod),iTPL+length(mod))
        patch(mod(iTPL).t_abs_Mod_pa,mod(iTPL).patch_norm_Mod,'g'); hold on
        plot(mod(iTPL).t_abs_Mod,mod(iTPL).y_norm_Mod,'k','linewidth',2)

        % Add Data Error Bars Normalized
        errorbar(mod(iTPL).t_abs_Mod,dat(iTPL).y_norm_Dat,...
            dat(iTPL).eps_norm_Dat(1,:),dat(iTPL).eps_norm_Dat(2,:),'linestyle','--','LineWidth',2);

        set(gca,'YScale','log','ylim',[0.01,3],'FontSize',16,'xtick',[0,20,40],'ytick',10.^[-2,-1,0])

    end


    %% Define Cell Counts at Each Time for Each Experiment
    % CellCounts = {[100,100,100,100],[100,100,100,100],[100,100,100,100,100],[100,100,100,100],[100,100,100,100]};
    dataFileDusp1 = 'Data/DUSP1_SSITcellresults_Final_Sep18.csv';
    allData = readtable(dataFileDusp1);
    allData = allData(allData.cyto_area_px>=12593&allData.cyto_area_px<=17685,:);
    allData = allData((allData.dex_conc == 100)|(allData.dex_conc == 0 & allData.time_tpl==0)|(allData.time == 0),:);
    for iTPL = 1:5
        for j = 1:length(mod(iTPL).t_abs_Mod)
            if j == 1
                allDataT = allData(isnan(allData.time_tpl)&allData.time==ttpl_list(iTPL),:);
            else
                allDataT = allData(allData.time_tpl==ttpl_list(iTPL)&...
                    allData.time==(ttpl_list(iTPL)+mod(iTPL).t_abs_Mod(j)),:);
            end
            CellCounts{iTPL}(j) = height(allDataT);
        end
    end


    %%
    % Pool Models and Data for Early and Late TPL Times
    EarlyMod.t_abs_Mod = mod(1).t_abs_Mod;
    EarlyMod.y_norm_Mod = (mod(1).y_norm_Mod.*CellCounts{1} + ...
        mod(2).y_norm_Mod.*CellCounts{2})./...
        (CellCounts{1}+CellCounts{2});
    EarlyMod.eps_norm_Mod = (mod(1).eps_norm_Mod.*[CellCounts{1},CellCounts{1}(end:-1:1)] + ...
        mod(2).eps_norm_Mod.*[CellCounts{2},CellCounts{2}(end:-1:1)])./...
        ([CellCounts{1},CellCounts{1}(end:-1:1)]+[CellCounts{2},CellCounts{2}(end:-1:1)]);
    EarlyMod.patch = [EarlyMod.y_norm_Mod,EarlyMod.y_norm_Mod(end:-1:1)]+...
        [-EarlyMod.eps_norm_Mod(1:end/2),EarlyMod.eps_norm_Mod(end/2+1:end)];
    EarlyMod.y_vals = [mod(1).y_norm_Mod(2:end), ...
        mod(2).y_norm_Mod(2:end)];
    
    
    EarlyDat.y_vals = [dat(1).y_norm_Dat(2:end), ...
        dat(2).y_norm_Dat(2:end)];
    EarlyDat.t_vals = [EarlyMod.t_abs_Mod(2:end),EarlyMod.t_abs_Mod(2:end)];

    EarlyDat.y_norm_Dat = (dat(1).y_norm_Dat.*CellCounts{1} + ...
        dat(2).y_norm_Dat.*CellCounts{2})./...
        (CellCounts{1}+CellCounts{2});
    EarlyDat.eps_norm_Dat = (dat(1).eps_norm_Dat.*CellCounts{1} + ...
        dat(2).eps_norm_Dat.*CellCounts{2})./...
        (CellCounts{1}+CellCounts{2});

    LateMod.t_abs_Mod = mod(3).t_abs_Mod;
    J = [1,3,4,5]; % The TPL at 150 and 180 min are missing the 10 minute time point
    LateMod.y_norm_Mod(J) = (mod(3).y_norm_Mod(J).*CellCounts{3}(J) + ...
        mod(4).y_norm_Mod.*CellCounts{4} + ...
        mod(5).y_norm_Mod.*CellCounts{5})./...
        (CellCounts{3}(J)+CellCounts{4}+CellCounts{5});
    LateMod.y_norm_Mod(2) = mod(3).y_norm_Mod(2);
    LateMod.eps_norm_Mod([1,3,4,5,6,7,8,10]) = (mod(3).eps_norm_Mod([1,3,4,5,6,7,8,10]).*[CellCounts{3}(J),CellCounts{3}(J(end:-1:1))] + ...
        mod(4).eps_norm_Mod.*[CellCounts{4},CellCounts{4}(end:-1:1)] + ...
        mod(5).eps_norm_Mod.*[CellCounts{5},CellCounts{5}(end:-1:1)])./...
        ([CellCounts{3}(J),CellCounts{3}(J(end:-1:1))]+[CellCounts{4},CellCounts{4}(end:-1:1)]+[CellCounts{5},CellCounts{5}(end:-1:1)]);
    LateMod.eps_norm_Mod([2,9]) = mod(3).eps_norm_Mod([2,9]);
    LateMod.patch = [LateMod.y_norm_Mod,LateMod.y_norm_Mod(end:-1:1)]+...
        [-LateMod.eps_norm_Mod(1:end/2),LateMod.eps_norm_Mod(end/2+1:end)];
    LateMod.y_vals = [mod(3).y_norm_Mod(2:end), ...
        mod(4).y_norm_Mod(2:end), mod(5).y_norm_Mod(2:end)];

    LateDat.y_vals = [dat(3).y_norm_Dat(2:end), ...
        dat(4).y_norm_Dat(2:end), dat(5).y_norm_Dat(2:end)];
    LateDat.t_vals = [LateMod.t_abs_Mod(2:end),EarlyMod.t_abs_Mod(2:end),EarlyMod.t_abs_Mod(2:end)];

    LateDat.y_norm_Dat(J) = (dat(3).y_norm_Dat(J).*CellCounts{3}(J) + ...
        dat(4).y_norm_Dat.*CellCounts{4} + ...
        dat(5).y_norm_Dat.*CellCounts{5})./...
        (CellCounts{3}(J)+CellCounts{4}+CellCounts{5});
    LateDat.y_norm_Dat(2) = dat(3).y_norm_Dat(2);
    LateDat.eps_norm_Dat(:,J) = (dat(3).eps_norm_Dat(:,J).*CellCounts{3}(J) + ...
        dat(4).eps_norm_Dat.*CellCounts{4} + ...
        dat(4).eps_norm_Dat.*CellCounts{5})./...
        (CellCounts{3}(J)+CellCounts{4}+CellCounts{5});
    LateDat.eps_norm_Dat(:,2) = dat(3).eps_norm_Dat(:,2);

    figure(101+10*kCompartment); clf;
    plot(EarlyMod.t_abs_Mod,EarlyMod.y_norm_Mod,'LineWidth',3); hold on;
    plot(EarlyMod.t_abs_Mod,EarlyDat.y_norm_Dat,'--ms','MarkerSize',16,'MarkerFaceColor','m'); hold on;
    plot(LateMod.t_abs_Mod,LateMod.y_norm_Mod,'LineWidth',3); hold on;
    plot(LateMod.t_abs_Mod,LateDat.y_norm_Dat,'--cs','MarkerSize',16,'MarkerFaceColor','c'); hold on;
    set(gca,'yscale','log','FontSize',16,'yLim',[0.08,2])

    f = figure(102+10*kCompartment); clf;
    set(f,'Name',[compartments{kCompartment},' DUSP1 vs Time after Triptolide (<=20 min Dex)']);
    patch([EarlyMod.t_abs_Mod,EarlyMod.t_abs_Mod(end:-1:1)],EarlyMod.patch,'g'); hold on;
    plot(EarlyMod.t_abs_Mod,EarlyMod.y_norm_Mod,'LineWidth',3); hold on;
    errorbar(EarlyMod.t_abs_Mod,EarlyDat.y_norm_Dat,EarlyDat.eps_norm_Dat(1,:),EarlyDat.eps_norm_Dat(2,:),'--ms','MarkerSize',16,'MarkerFaceColor','m'); hold on;
    set(gca,'yscale','log','FontSize',16,'yLim',[0.09,3])

    f = figure(103+10*kCompartment); clf;
    set(f,'Name',[compartments{kCompartment},' DUSP1 vs Time after Triptolide (>=75 min Dex)']);
    patch([LateMod.t_abs_Mod,LateMod.t_abs_Mod(end:-1:1)],LateMod.patch,'g'); hold on;
    plot(LateMod.t_abs_Mod,LateMod.y_norm_Mod,'LineWidth',3); hold on;
    errorbar(LateMod.t_abs_Mod,LateDat.y_norm_Dat,LateDat.eps_norm_Dat(1,:),LateDat.eps_norm_Dat(2,:),'--ms','MarkerSize',16,'MarkerFaceColor','m'); hold on;
    set(gca,'yscale','log','FontSize',16,'yLim',[0.09,3])

end
%%
[tstat, p, mEarly, mLate, semEarly, semLate] = compare_slopes(EarlyDat.t_vals,log(EarlyDat.y_vals),LateDat.t_vals,log(LateDat.y_vals),'right');
disp('DATA')
disp(['Early Decay Slope: ',num2str(mEarly),' +/- ',num2str(semEarly)])

% [mLate,semLate,bLate] = slope_sem(LateDat.t_vals,log(LateDat.y_vals));
disp(['Late Decay Slope: ',num2str(mLate),' +/- ',num2str(semLate)])

disp(['1-sided T-test, t = ', num2str(tstat)])
disp(['                p = ', num2str(p)])

disp('MODEL')
[tstat, p, mEarly, mLate, semEarly, semLate] = compare_slopes(EarlyDat.t_vals,log(EarlyMod.y_vals),LateDat.t_vals,log(LateMod.y_vals),'right');
disp(['Early Decay Slope: ',num2str(mEarly),' +/- ',num2str(semEarly)])

% [mLate,semLate,bLate] = slope_sem(LateDat.t_vals,log(LateDat.y_vals));
disp(['Late Decay Slope: ',num2str(mLate),' +/- ',num2str(semLate)])

disp(['1-sided T-test, t = ', num2str(tstat)])
disp(['                p = ', num2str(p)])


function [tstat, p, m1, m2, sem1, sem2] = compare_slopes(x1,y1,x2,y2,tail)
% Compare slopes of two regression lines using t-test
%
% Outputs:
%   tstat - t statistic
%   p     - two-sided p-value
%   m1,m2 - slopes
%   sem1,sem2 - standard errors of slopes

    % --- dataset 1 ---
    x1 = x1(:);
    y1 = y1(:);

    valid1 = ~(isnan(x1) | isnan(y1));
    x1 = x1(valid1);
    y1 = y1(valid1);

    n1 = length(x1);

    pfit1 = polyfit(x1,y1,1);
    m1 = pfit1(1);

    yfit1 = polyval(pfit1,x1);
    resid1 = y1 - yfit1;

    s2_1 = sum(resid1.^2)/(n1-2);
    sem1 = sqrt(s2_1 / sum((x1-mean(x1)).^2));

    % --- dataset 2 ---
    x2 = x2(:);
    y2 = y2(:);

    valid2 = ~(isnan(x2) | isnan(y2));
    x2 = x2(valid2);
    y2 = y2(valid2);

    n2 = length(x2);

    pfit2 = polyfit(x2,y2,1);
    m2 = pfit2(1);

    yfit2 = polyval(pfit2,x2);
    resid2 = y2 - yfit2;

    s2_2 = sum(resid2.^2)/(n2-2);
    sem2 = sqrt(s2_2 / sum((x2-mean(x2)).^2));

    % --- compare slopes ---
    tstat = (m1 - m2) / sqrt(sem1^2 + sem2^2);

    df = n1 + n2 - 4;

    % One-sided p-value
    switch lower(tail)

        case 'right'   % H1: m1 > m2
            p = 1 - tcdf(tstat, df);

        case 'left'    % H1: m1 < m2
            p = tcdf(tstat, df);

        otherwise
            error('tail must be ''right'' or ''left''');
    end
end