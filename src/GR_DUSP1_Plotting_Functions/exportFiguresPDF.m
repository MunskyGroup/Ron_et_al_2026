close all
files = dir('savedFigures/FinalFigures/MainFigures/');
files = files(contains({files.name},'.fig'));

for i = 1:length(files)
    fig = openfig(['savedFigures/FinalFigures/MainFigures/',files(i).name]);
    exportgraphics(fig,['savedFigures/FinalFigures/MainFigures/pdf/',files(i).name(1:end-3),'pdf'])
end

%%
close all
files = dir('savedFigures/FinalFigures/SupplementalFigures/');
files = files(contains({files.name},'.fig'));

for i = 1:length(files)
    fig = openfig(['savedFigures/FinalFigures/SupplementalFigures/',files(i).name]);
    exportgraphics(fig,['savedFigures/FinalFigures/SupplementalFigures/pdf/',files(i).name(1:end-3),'pdf'])
end