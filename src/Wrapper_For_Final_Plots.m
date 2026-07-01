function [Models] = Wrapper_For_Final_Plots(opts)
%% Generate fits for GR and nuclear DUSP1 at 100nM Dex
% Quantify all of the fits for DEX and titrations, with optional arguments:
%   * refitGR = false
%   * refitDUSP1nuc = false
%   * refitDUSP1cyt = false
%   * refitKTPTL = false
%   * whichGRmodel = []
%   * whichDUSP1modelnuc = []
%   * whichDUSP1modelcyt = []
%   * makePlots = false
%   * finalDUSP1model = []
%   * makePredictions = false
%   * useCluster = false
%   * recomputeSSA = false

% GR Models 1-4:
%   * Base Model: Nuc + Cyt Degradation, Dex -> Cyt2Nuc
%   * Alt Model:  Cyt Degradation, Dex -> Cyt2Nuc
%   * Alt Model:  Nuc Degradation, Dex -> Cyt2Nuc
%   * Alt Model:  Nuc + Cyt Degradation, Dex -| Nuc2Cyt

% DUSP1 Nuclear Models 1-8:
%   * case 1 % KOFF
%         -- Free Parameters: 
%               koff, kon_0, kr_on_0, knuc2cyt, tau_elong, mkoff
%         -- Sets to zero: 
%               kon_1, kr_on_1
%   * case 2 % KON
%         -- Free Parameters: 
%               koff, kon_0, kr_on_0, knuc2cyt, tau_elong, kon_1
%         -- Sets to zero: 
%               mkoff, kr_on_1
%   * case 3 % KR
%         -- Free Parameters: 
%               koff, kon_0, kr_on_0, knuc2cyt, tau_elong, kr_on_1
%         -- Sets to zero: 
%               kon_1, mkoff
%   * case 4 % General 2-State (KON, KOFF, KR)
%         -- Free Parameters: 
%               koff, kon_0, kr_on_0, knuc2cyt, tau_elong, kr_off, kon_1,
%               mkoff, kr_on_1
%         -- Sets to zero: 
%               
%   * case 5 % General 2-State (KON, KOFF)
%         -- Free Parameters: 
%               koff, kon_0, kr_on_0, knuc2cyt, tau_elong, kr_off, kon_1,
%               mkoff
%         -- Sets to zero: 
%               kr_on_1
%   * case 6 % General 2-State (KON, KOFF, KR) - no production in OFF state
%         -- Free Parameters: 
%               koff, kon_0, kr_on_0, knuc2cyt, tau_elong, kon_1, mkoff,
%               kr_on_1
%         -- Sets to zero: 
%               kr_off
%   * case 7 % General 2-State (KON, KOFF) - no production in OFF state
%         -- Free Parameters: 
%               koff, kon_0, kr_on_0, knuc2cyt, tau_elong, kon_1, mkoff
%         -- Sets to zero: 
%               kr_off, kr_on_1
%   * case 8 % General 2-State (KON) - no production in OFF state
%         -- Free Parameters: 
%               koff, kon_0, kr_on_0, knuc2cyt, tau_elong, kon_1 
%         -- Sets to zero: 
%               kr_off, mkoff, kr_on_1
% DUSP1 Cytoplasmic ("Final") Models 1-12, except 8-10 (deprecated models):
%   * case 1 % Kon, Koff, Kr free, 
%              use Nuc fit parameters and fit only Cyt Parameters
%   *     -- Free Parameters:
%               degCyt0, degCyt1, kDTTP, TTP1, tTTP, eta
%   *     -- Sets to zero:
%               kr_off
%   *     -- Does not include TS (Cyt parameters do not affect the TS)
%   * case 2 % Kon, Koff, Kr free
%   *     -- Free Parameters:
%               koff, kon_0, kr_on_0, knuc2cyt, degCyt0, degCyt1, kDTTP, 
%               tau_long, TTP1, tTTP, eta, kon_1, mkoff, kr_on_1
%   *     -- Sets to zero:
%               kr_off
%   * case 3 % No Saturation
%   *     -- Free Parameters:
%               koff, kon_0, kr_on_0, knuc2cyt, degCyt0, degCyt1, 
%               tau_long, TTP1, tTTP, eta, kon_1, mkoff, kr_on_1
%   *     -- Sets to zero:
%               kr_off
%   * case 4 % No Time Variation
%   *     indsFreePars =  [1,2,3,4,14,15,16,19,20,21,24,25,26];
%   *     -- Free Parameters:
%               koff, kon_0, kr_on_0, knuc2cyt, degCyt0, degCyt1, kDTTP, 
%               tau_long, TTP1, tTTP, kon_1, mkoff, kr_on_1
%   *     -- Sets to zero:
%               kr_off
%   * case 5 % KON Only
%   *     -- Free Parameters:
%               koff, kon_0, kr_on_0, knuc2cyt, degCyt0, degCyt1, kDTTP, 
%               tau_long, TTP1, tTTP, eta, kon_1
%   *     -- Sets to zero:
%               kr_off, mkoff, kr_on_1
%   * case 6 % KON Only - No Prior (same as case 5 but without a prior)
%   *     indsFreePars =  [1,2,3,4,14,15,16,19,20,21,22,24];
%   *     -- Free Parameters:
%               koff, kon_0, kr_on_0, knuc2cyt, degCyt0, degCyt1, kDTTP, 
%               tau_long, TTP1, tTTP, eta, kon_1
%   *     -- Sets to zero:
%               kr_off, mkoff, kr_on_1
%   * case 7 % Full - No Prior (same as case 2 but without a prior)
%   *     -- Free Parameters:
%               koff, kon_0, kr_on_0, knuc2cyt, degCyt0, degCyt1, kDTTP, 
%               tau_long, TTP1, tTTP, eta, kon_1, mkoff, kr_on_1
%   *     -- Sets to zero:
%               kr_off
%   * case 8 % DO NOT USE
%   * case 9 % DO NOT USE    
%   * case 10 % DO NOT USE
%   * case 11  %% Mechanistic model with TTP
%   *     -- Free Parameters:
%               koff, kon_0, kr_on_0, knuc2cyt, degCyt0, degCyt1,  
%               tau_long, TTP1, tTTP, eta, kon_1, mkoff, kr_on_1,
%               kon_1_TTP, koff_TTP, kr_on_TTP, knuc2cyt_TTP, ktranslate, 
%               degProt, kbind, TTP_eta, unbind, degCyt0_TTP
%   *     -- Sets to zero:
%               kr_off
%   * case 12  %% Mechanistic model with TTP and P38
%   *     -- Free Parameters:
%               koff, kon_0, kr_on_0, knuc2cyt, degCyt0, degCyt1,  
%               tau_long, TTP1, tTTP, eta, kon_1, mkoff, kr_on_1,
%               kon_1_TTP, koff_TTP, kr_on_TTP, knuc2cyt_TTP, ktranslate, 
%               degProt, kbind, TTP_eta, unbind, degCyt0_TTP,     
%               kMPK1, gMPK1, kP38, gP38, pP38, dP38, sequest, release
%   *     -- Sets to zero:
%               kr_off

makePlots = isfield(opts,'remakeFigs')&&opts.remakeFigs;
refitMods = isfield(opts,'refitModels')&&opts.refitModels;

if isfield(opts,'regenLibrary')&&opts.regenLibrary
    disp('********************************************************')
    disp('Deleting Stored Model Library')
    disp('********************************************************')
    delete savedParameters/GRDusp1ModelLibrary.mat
end

if isfield(opts,'hideFigures')&&opts.hideFigures
    set(0,'defaultFigureVisible','off')
    hideFigures = true;
else
    set(0,'defaultFigureVisible','on')
    hideFigures = false;
end

if ~isfield(opts,'recomputeSSA')
    opts.recomputeSSA = false;
end
recomputeSSA = opts.recomputeSSA;


%% Make GR Plots
GRMods = [];
if isfield(opts,'GRModels')
    choices = {'DexControlsGRImport_NucCytDeg',...
        'DexControlsGRImport_NucDeg','DexControlsGRImport_CytDeg',...
        'DexControlsGRExport_NucCytDeg'};
    if strcmpi(opts.GRModels{1},'all')
        GRMods = [1:4];
    else
        GRMods = find(contains(choices,opts.GRModels));
    end    
end
if ~isempty(GRMods)
    disp('********************************************************')
    disp('Running Analyses for GR Models')
    disp('********************************************************')
    GR_DUSP1_Fitting_Final(whichGRmodel=GRMods, makePlots=makePlots, ...
        refitGR=refitMods, recomputeGR=true);
    close all
    for iMod = 1:length(GRMods)
        disp('********************************************************')
        disp(['Making Plots for GR Model: ',choices{iMod}])
        disp('********************************************************')
        makeFinalFigures(plots = {'GRdistributions','GRModelComparison', 'TotalGR'}, ...
            hideFigures=hideFigures, GRmodel=iMod);
        close all
    end
end


%% Make DUSP1 Nuc Fitting Plots
% SKIPPED
% GR_DUSP1_Fitting_Final(whichDUSP1modelnuc=1:8, makePlots=true)

%% Make Nuc/Cyt Fitting Plots
% GR_DUSP1_Fitting_Final(whichDUSP1modelcyt=[3,4,7,8,11], ...
%                         recomputeSSA=true, makePlots=true)

%% Generate fits and predictions for DUSP1Nuc/Cyt at all Dex levels
% Parameters = GR_DUSP1_Fitting_Final(whichDUSP1modelcyt=[3,4,7,8,11],...
%                                       makePredictions=true, ...
%                                       recomputeSSA=true)

%% Make and Save Prediction Plots for Final Models
FullMods = [];
if isfield(opts,'FullModels')
    if strcmpi(opts.FullModels{1},'all')
        FullMods = [1:5];
    else
        choices = {'TVandSat','NoTV','NoSat','LinDecay','Mechanistic'};
        FullMods = find(contains(choices,opts.FullModels));
    end    
end
if ~isempty(FullMods)    
    for i = 1:length(FullMods)
        switch FullMods(i)
            case 1
                disp('********************************************************')
                disp('Running Analyses for Full Final Semi-Mechanistic Model')
                disp('********************************************************')
                finalDUSP1model=2;
           case 2
                disp('********************************************************')
                disp('Running Analyses for Full No Time Variation Model')
                disp('********************************************************')
                finalDUSP1model=4;
           case 3
                disp('********************************************************')
                disp('Running Analyses for Full No Saturation Model')
                disp('********************************************************')
                finalDUSP1model=3;
           case 4
                disp('********************************************************')
                disp('Running Analyses for Full First Order Decay Model')
                disp('********************************************************')
                finalDUSP1model=8;
           case 5
                disp('********************************************************')
                disp('Running Analyses for Full Extended Mechanistic Model')
                disp('********************************************************')
                finalDUSP1model=10;
        end
        close all
        GR_DUSP1_Fitting_Final(refitDUSP1cyt=refitMods, ...
            finalDUSP1model=finalDUSP1model, whichDUSP1modelcyt=finalDUSP1model , makePlots=makePlots, ...
            recomputeSSA=recomputeSSA, makeTPLPlots=makePlots, ...
            makeTSPlots=makePlots, makeTitrPlots=makePlots, ...
            hidefigures=hideFigures, makePredictions=true);
        close all

        makeFinalFigures(finalModel=finalDUSP1model, hideFigures=hideFigures);

    end
end

disp('********************************************************')
disp('Generating MH Parameter Uncertainty Tables')
disp('********************************************************')
makeFinalFigures(plots = {'MHGR','MHCyt','MHMech'});