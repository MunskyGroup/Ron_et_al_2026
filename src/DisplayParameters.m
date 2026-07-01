function [] = DisplayParameters()
modelLibrary = 'savedParameters/GRDusp1ModelLibrary';
dataFileGR = 'Data_Aug2025/GR_SSITcellresults_Final_Sep18.csv';
dataFileDusp1 = 'Data_Aug2025/DUSP1_SSITcellresults_Final_Sep18.csv';
savedParsFile = 'savedParameters/Parameters_Final';

% Full display index list: includes GR shared params (5:12) and kDTTP (16)
displayInds = [1,2,3,4, 5,6,7,8,9,10,11,12, 14,15,19,20,21,22,24,25,26,16,27:36,18];
fixedInds = 5:12;  % GR pars: fixed from GR fit, never NaN

%% Load GR Parameters
combined_GRModel = dusp1ModelLibrary_Final('combinedGRModel',false,modelLibrary,dataFileGR,dataFileDusp1);
load(savedParsFile,'GRpars');
combined_GRModel = combined_GRModel.updateModels(GRpars, false);
GRparValues = combined_GRModel.SSITModels{3}.parameters(5:12, 2);

%% DUSP1 TS/Nuc Parameters -- KON only model
disp('TS/Nuc Model - KON only, no expression in off state')
DUSP1_Model_100nM = dusp1ModelLibrary_Final('ModelDUSP1_100nM',false,modelLibrary);
DUSP1_Model_100nM.parameters(5:12, 2) = GRparValues;
load(savedParsFile,'Pars_KON_zeroOFF');
modelVarsToFitA = [1,2,3,4,19,24];
namInds = setdiff(displayInds, [modelVarsToFitA, fixedInds]);
DUSP1_Model_100nM.parameters(modelVarsToFitA, 2) = num2cell(Pars_KON_zeroOFF);
DUSP1_Model_100nM.parameters(namInds, 2) = {NaN};
T.NucKON = DUSP1_Model_100nM.parameters(displayInds, 2);

%% DUSP1 TS/Nuc Parameters
disp('TS/Nuc Model - all parameters free, no expression in off state')
DUSP1_Model_100nM = dusp1ModelLibrary_Final('ModelDUSP1_100nM',false,modelLibrary);
DUSP1_Model_100nM.parameters(5:12, 2) = GRparValues;
load(savedParsFile,'Pars_KON_KOFF_KR_zeroOFF');
modelVarsToFitA = [1,2,3,4,19,24,25,26];
namInds = setdiff(displayInds, [modelVarsToFitA, fixedInds]);
DUSP1_Model_100nM.parameters(modelVarsToFitA, 2) = num2cell(Pars_KON_KOFF_KR_zeroOFF);
DUSP1_Model_100nM.parameters(namInds, 2) = {NaN};
T.NucALL = DUSP1_Model_100nM.parameters(displayInds, 2);

%% DUSP1 TS/Nuc/Cyt Parameters
disp('TS/Nuc/Cyt Model - all parameters free, no expression in off state')
SSAModel_100 = dusp1ModelLibrary_Final('fullSSAModel_100_2StateGen',false,modelLibrary);
SSAModel_100.parameters(5:12, 2) = GRparValues;
load(savedParsFile,'FinalModel_allButGR');
load(savedParsFile,'ktptl_FinalModel_allButGR')
indsFreePars_NucCyt = [1,2,3,4,14,15,16,19,20,21,22,24,25,26];
namInds = setdiff(displayInds, [indsFreePars_NucCyt, fixedInds]);
SSAModel_100.parameters(indsFreePars_NucCyt, 2) = num2cell(FinalModel_allButGR);
SSAModel_100.parameters(namInds, 2) = {NaN};
SSAModel_100.parameters(18, 2) = {ktptl_FinalModel_allButGR};
T.NucCyt = SSAModel_100.parameters(displayInds, 2);

% disp('TS/Nuc/Cyt Model - all parameters free, no repression in off state')
% SSAModel_100 = dusp1ModelLibrary_Final('fullSSAModel_100_2StateGen',false,modelLibrary);
% SSAModel_100.parameters(5:12, 2) = GRparValues;
% load(savedParsFile,'KONModel_allButGR');
% load(savedParsFile,'ktptl_KONModel_allButGR')
% indsFreePars_NucCyt = [1,2,3,4,14,15,16,19,20,21,22,24];
% namInds = setdiff(displayInds, [indsFreePars_NucCyt, fixedInds]);
% SSAModel_100.parameters(indsFreePars_NucCyt, 2) = num2cell(KONModel_allButGR);
% SSAModel_100.parameters(namInds, 2) = {NaN};
% SSAModel_100.parameters(18, 2) = {ktptl_KONModel_allButGR};
% T.NucCyt_KonOnly = SSAModel_100.parameters(displayInds, 2);

%% Extended TTP Model
% Note: ITTP expression is removed in this model (propensity 10 = 'degCyt0*rCyt')
% TTP1 (20), tTTP (21), and eta (22) are structurally unused but were included
% during fitting.
% kDTTP (16) is NaN as it only appears in ITTP.
disp('Extended Model with TTP')
SSAModel_TTP = dusp1ModelLibrary_Final('GR_DUSP1_TTP_Model_2',false,modelLibrary);
SSAModel_TTP.parameters(5:12, 2) = GRparValues;
load(savedParsFile,'Pars_GR_DUSP1_TTP_Model_All');
load(savedParsFile,'ktptl_GR_DUSP1_TTP_Model_All')
indsFreePars_TTP = [1,2,3,4,14,15,19,24,25,26,27:36];  % 20 elements -- matches saved pars
assert(numel(Pars_GR_DUSP1_TTP_Model_All)==numel(indsFreePars_TTP), ...
    'Mismatch: Pars_GR_DUSP1_TTP_Model_All has %d elements but indsFreePars_TTP has %d', ...
    numel(Pars_GR_DUSP1_TTP_Model_All), numel(indsFreePars_TTP))
namInds = [setdiff(displayInds, [indsFreePars_TTP, fixedInds])];
SSAModel_TTP.parameters(indsFreePars_TTP, 2) = num2cell(Pars_GR_DUSP1_TTP_Model_All);
SSAModel_TTP.parameters(namInds, 2) = {NaN};
SSAModel_TTP.parameters(18, 2) = {ktptl_GR_DUSP1_TTP_Model_All};
T.Extended = SSAModel_TTP.parameters(displayInds, 2);

%% Diagnostics
disp('=== Parameter Assignment Diagnostics ===')

% Helper to check expected NaN/non-NaN assignments
function checkPars(modelName, params, shouldBeNaN, shouldBeSet)
    fprintf('\n-- %s --\n', modelName)
    for i = 1:numel(shouldBeNaN)
        idx = shouldBeNaN(i);
        val = params{idx,2};
        if isnan(val); status = 'OK'; else; status = 'FAIL (expected NaN)'; end
        fprintf('  [%s] %s (idx %d): %s\n', status, params{idx,1}, idx, num2str(val))
    end
    for i = 1:numel(shouldBeSet)
        idx = shouldBeSet(i);
        val = params{idx,2};
        if ~isnan(val); status = 'OK'; else; status = 'FAIL (expected value)'; end
        fprintf('  [%s] %s (idx %d): %s\n', status, params{idx,1}, idx, num2str(val))
    end
end

% TS/Nuc: cyto deg, kDTTP, TTP1, tTTP, eta should be NaN; fitted + GR fixed should be set
checkPars('TS/Nuc', DUSP1_Model_100nM.parameters, ...
    [14,15,16,20,21,22], ...           % expected NaN: degCyt0/1, kDTTP, TTP1, tTTP, eta
    [modelVarsToFitA, fixedInds])      % expected set: fitted + GR

% TS/Nuc/Cyt: extended TTP species params (27:36) should be NaN; fitted + GR fixed should be set
checkPars('TS/Nuc/Cyt', SSAModel_100.parameters, ...
    [27:36], ...                       % expected NaN: extended TTP model params
    [indsFreePars_NucCyt, fixedInds])  % expected set: fitted + GR

% Extended: kDTTP (16) should be NaN; TTP1/tTTP/eta present but structurally unused
checkPars('Extended TTP', SSAModel_TTP.parameters, ...
    [16], ...                          % expected NaN: kDTTP only
    [indsFreePars_TTP, fixedInds])     % expected set: fitted + GR

% Verify TTP_eta exists and is correctly assigned in extended model
ttpEtaIdx = find(strcmp(SSAModel_TTP.parameters(:,1), 'TTP_eta'));
if ~isempty(ttpEtaIdx)
    fprintf('\nTTP_eta found at index %d, value = %s\n', ttpEtaIdx, ...
        num2str(SSAModel_TTP.parameters{ttpEtaIdx,2}))
    if ~ismember(ttpEtaIdx, indsFreePars_TTP)
        warning('TTP_eta (idx %d) is NOT in indsFreePars_TTP — value was not updated from saved parameters', ttpEtaIdx)
    end
else
    warning('TTP_eta not found in extended model parameter list')
end

%% Assemble Table
T.Names = SSAModel_TTP.parameters(displayInds, 1);

% Tab = table(T.Names, T.Nuc, T.NucCyt, T.NucCyt_KonOnly, T.Extended);
Tab = table(T.Names, T.NucALL, T.NucKON, T.NucCyt, T.Extended);
Tab.Properties.VariableNames = {'Parameters','TS/Nuc (All free)','TS/Nuc (KON Only)','TS/Nuc/Cyt','Extended'};

%% Reorder Table
order = [1,2,19,20,3,21,15,4,13,14,16,17,18,22,23,24,25,26,32,27,28,29,31,30];
orderedTable = Tab(order,:);

%% Write to LaTeX
fid = fopen('ParameterTables/ParameterComparisons.txt','w');
fprintf(fid,'\n%%Comparison of model parameter values:\n');
fprintf(fid,'\\begin{tabular}{|l|l|l|l|l|l|}\n');
fprintf(fid,'\\hline\n');
fprintf(fid,'\\textbf{Description (units)} & \\textbf{Parameter} & \\makecell{\\textbf{TS/Nuc}\\\\ ($k_{on},k_{off},k_{r}$)}& \\makecell{\\textbf{TS/Nuc} \\\\($k_{on}$ only)}& \\textbf{TS/Nuc/Cyt} & \\textbf{Extended} \\\\ \\hline\n');

% fprintf(fid,'Parameter & TS/Nuc ($k_{on}$,$k_{off}$,$k_{r}$) & TS/Nuc ($k_{on}$ only) & TS/Nuc/Cyt ($k_{on}$,$k_{off}$,$k_{r}$) & Mechanistic Model \\hline \\\\ \n');

descriptions = {'Basal gene inactivation rate (min$^{-1}$)' ,'$k_{off}$'        ;
    'Basal gene activation rate (min$^{-1}$)'              ,'$k_{on,0}$';
    'nucGR-dependent gene activation rate (min$^{-1}$)'    ,'$k_{on,1}$';
    'nucGR modulation of inactivation (molecules$^{-1}$)'  ,'$m_{koff}$';
    'Basal DUSP1 transcription rate (min$^{-1}$)'          ,'$kr_{on,0}$';
    'nucGR-dependent DUSP1 transcription rate (min$^{-1}$)','$kr_{on,1}$';
    'DUSP1 elongation time (min)'                          ,'$\tau_{elong}$';
    'DUSP1 nuclear export rate (min$^{-1}$)'               ,'$k_{nc}$';
    'Basal DUSP1 mRNA degradation rate (min$^{-1}$)'       ,'$\gamma_{Cyt,0}$';
    'Max TTP-dependent degradation rate (min$^{-1}$)'      ,'$\gamma_{Cyt,1}$';
    'Max fold-change in TTP-mediated degradation'          ,'$TTP1$';
    'Half-max time for TTP activity (min)'                 ,'$t_{TTP}$';
    'Hill coefficient, TTP time-dependence'                ,'$\eta$';
    'TTP-mRNA binding saturation (molecules)'              ,'$kD_{TTP}$';
    'nucGR-dependent TTP gene activation rate (min$^{-1}$)','$kon_{TTP,1}$';
    'TTP gene inactivation rate (min$^{-1}$)'              ,'$koff_{TTP}$';
    'TTP transcription rate (min$^{-1}$)'                  ,'$kr_{on,TTP}$';
    'TTP nuclear export rate (min$^{-1}$)'                 ,'$k_{TTP,nc}$';
    'Basal TTP mRNA degradation rate (min$^{-1}$)'         ,'$\gamma_{TTP,Cyt}$';
    'TTP translation rate (min$^{-1}$)'                    ,'$k_{translate}$';
    'TTP protein degradation rate (min$^{-1}$)'            ,'$\gamma_{CTT,Prot}$';
    'TTP–mRNA binding rate (min$^{-1}$ mol$^{-1}$)'        ,'$k_{bind}$';
    'TTP–mRNA unbinding rate (min$^{-1}$)'                 ,'$k_{unbind}$';
    'Hill coefficient, TTP–mRNA binding'                   ,'$\eta_{TTP}$'     };

for p = 1:size(orderedTable,1)

    vals = cellfun(@(x) ...
        ternary(isnan(x), '-', sprintf('%.2e', x)), ...
        orderedTable{p,2:5}, ...
        'UniformOutput', false);

    fprintf(fid,'%15s & %15s & %s & %s & %s & %s \\\\ \\hline\n', ...
       descriptions{p,1}, descriptions{p,2}, vals{:});

end
fprintf(fid,'\\end{tabular}');
fclose(fid);

end

function out = ternary(cond,a,b)
    if cond
        out = a;
    else
        out = b;
    end
end


