function figHandles = makeGRPlots(combinedModel,GRpars,opts)
arguments
    combinedModel
    GRpars
    opts.splitReps = false
    opts.GR_Data = 'Data/GR_SSITcellresults_Final_Sep18.csv'
end
splitReps = opts.splitReps;
GR_Data = opts.GR_Data;

GRfitCases = {'1','1',101,'GR Fit (1nM Dex)';...
    '10','10',102,'GR Fit (10nM Dex)';...
    '100','100',103,'GR Fit (100nM Dex)'};

combinedGRModel = combinedModel.updateModels(GRpars,false);
nMods = length(combinedGRModel.SSITModels);
ModelGroup = cell(nMods,1);
figHandles = cell(nMods, 1); % Initialize figHandles to store figure handles
for i=1:nMods
    if splitReps
        for rep = {'A','B','C'}
            %  Update parameters in original models.
            ModelGroup{i} = combinedGRModel.SSITModels{i};
            ModelGroup{i} = ModelGroup{i}.loadData(GR_Data,...
                {'nucGR','normGRnuc';'cytGR','normGRcyt'},...
                {[],[], ...
                ['(TAB.dex_conc==',GRfitCases{i,1},...
                '|(TAB.dex_conc==0&TAB.time==0&TAB.imageDates~=20230522))&TAB.time~=20' ...
                '&TAB.time~=40&TAB.time~=60' ...
                '&TAB.time~=90&TAB.time~=150' ...
                '&strcmp(TAB.replica,''',rep{1},''')']});
            ModelGroup{i}.tSpan = sort(unique([ModelGroup{i}.tSpan,linspace(0,180,30)]));
            if strcmp(rep,'A')
                figHandles{i} = ModelGroup{i}.makeFitPlot([],1,[],true,'IQR',0.25,true);
            else
                figHandles{i} = ModelGroup{i}.makeFitPlot([],1,figHandles{i},true,'IQR',0.25,true);
            end
        end
    else
        %  Update parameters in original models.
        ModelGroup{i} = combinedGRModel.SSITModels{i};
        ModelGroup{i} = ModelGroup{i}.loadData(GR_Data,...
            {'nucGR','normGRnuc';'cytGR','normGRcyt'},...
            {[],[], ...
            ['(TAB.dex_conc==',GRfitCases{i,1},...
            '|(TAB.dex_conc==0&TAB.time==0&TAB.imageDates~=20230522))&TAB.time~=20' ...
            '&TAB.time~=40&TAB.time~=60' ...
            '&TAB.time~=90&TAB.time~=150']});
        ModelGroup{i}.tSpan = sort(unique([ModelGroup{i}.tSpan,linspace(0,180,30)]));
        figHandles{i} = ModelGroup{i}.makeFitPlot([],1,[],true,'IQR',0.25,false);
    end
end
end