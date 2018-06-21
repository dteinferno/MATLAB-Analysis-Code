clear;
clc;
ROIs = [1:16];
pkDiff = 3;
tileAmp = [0.1:0.1:1];
wedgeAmp = ones(1,length(tileAmp));
wedgeOffset = 7;
wedgeAct = cos(2*pi*(ROIs+wedgeOffset)/length(ROIs))+1;
% wedgeAct(find(wedgeAct<0)) = 0;
PBtileActL = cos(2*pi*(ROIs+wedgeOffset+pkDiff)/length(ROIs))+1;
PBtileActR = cos(2*pi*(ROIs+wedgeOffset-pkDiff)/length(ROIs))+1;
% PBtileAct(find(wedgeAct<0)) = 0;
EBtileActL = tileAmp'*PBtileActL+wedgeAmp'*wedgeAct;
EBtileActL = EBtileActL - 2;
EBtileActL(find(EBtileActL<0)) = 0;
EBtileActR = tileAmp'*PBtileActR+wedgeAmp'*wedgeAct;
EBtileActR = EBtileActR - 2;
EBtileActR(find(EBtileActR<0)) = 0;

figure;
hold on;
plot(EBtileActL','r')
plot(EBtileActR','r')
plot(wedgeAct','g','LineWidth',2)
xlim([1 16]);

