close all
A = load('savedResults/SSALikelihood/SeqABC_TVandSat_chain.mat');
B = load('savedResults/SSALikelihood/SeqABC_NoTV_chain.mat');
C = load('savedResults/SSALikelihood/SeqABC_NoSat_chain.mat');
D = load('savedResults/SSALikelihood/SeqABC_LinDecay_chain.mat');
E = load('savedResults/SSALikelihood/SeqABC_Mechanistic_chain.mat');
figure
for iii = 1:6
    for jjj = iii+1:7
        subplot(6,6,(iii-1)*6+jjj-1)
        scatter(A.smpl(~isnan(A.value),jjj),smpl(~isnan(A.value),iii),16,[1:sum(~isnan(A.value))])
        hold on
        scatter(B.smpl(~isnan(B.value),jjj),smpl(~isnan(B.value),iii),16,[1:sum(~isnan(B.value))])
        scatter(C.smpl(~isnan(C.value),jjj),smpl(~isnan(C.value),iii),16,[1:sum(~isnan(C.value))])
        scatter(D.smpl(~isnan(D.value),jjj),smpl(~isnan(D.value),iii),16,[1:sum(~isnan(D.value))])
        scatter(E.smpl(~isnan(E.value),jjj),smpl(~isnan(E.value),iii),16,[1:sum(~isnan(E.value))])
    end
end
figure
plot(A.value(~isnan(A.value)))
hold on
plot(B.value(~isnan(B.value)))
plot(C.value(~isnan(C.value)))
plot(D.value(~isnan(D.value)))
plot(E.value(~isnan(E.value)))
legend
