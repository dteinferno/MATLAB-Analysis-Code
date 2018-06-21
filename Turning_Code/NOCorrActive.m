GtOn = 0.42/log(2);
GtOff = 0.33/log(2);

CCs = {};

tStop = 5;
vRThresh = 0.000075;
numPts = 17;

% Plot overviews for each dataset
for flyID = 1:length(allFlyData)
    
    % Plot parameters for the dark data
    for darkID = 1:length(allFlyData{flyID}.dark);
        ROISig = allFlyData{flyID}.dark{darkID}.GROIave(1,:)-allFlyData{flyID}.dark{darkID}.GROIave(2,:);
        
        vRot = allFlyData{flyID}.dark{darkID}.positionDatMatch.vRot(1:end);
        tAll = allFlyData{flyID}.dark{darkID}.positionDatMatch.OffsetRotMatch(:,1);
        
        vRotCCW = vRot;
        vRotCCW(find(vRotCCW<0)) = 0;
        vRotCW = -vRot;
        vRotCW(find(vRotCW<0)) = 0;
        
        vRotCCWGoingUp = find(diff(vRotCCW)>0);
        vRotCCWGoingDown = find(diff(vRotCCW)<0);
        vRotCCWUp = zeros(length(vRotCCW),1);
        vRotCCWUp(vRotCCWGoingUp)= vRotCCW(vRotCCWGoingUp);
        vRotCCWDown = zeros(length(vRotCCW),1);
        vRotCCWDown(vRotCCWGoingDown)= vRotCCW(vRotCCWGoingDown);
        vRotCCWConv = conv(vRotCCWUp,exp(-tAll/GtOn)) + conv(vRotCCWDown,exp(-tAll/GtOff));
        
        vRotCWGoingUp = find(diff(vRotCW)>0);
        vRotCWGoingDown = find(diff(vRotCW)<0);
        vRotCWUp = zeros(length(vRotCW),1);
        vRotCWUp(vRotCWGoingUp)= vRotCW(vRotCWGoingUp);
        vRotCWDown = zeros(length(vRotCW),1);
        vRotCWDown(vRotCWGoingDown)= vRotCW(vRotCWGoingDown);
        vRotCWConv = conv(vRotCWUp,exp(-tAll/GtOn)) + conv(vRotCWDown,exp(-tAll/GtOff));
        
        vRotConv = vRotCCWConv-vRotCWConv;
        
        %find periods where the rotation is above some threshold
        vZero = find(abs(vRotConv(1:length(vRot)))<vRThresh);
        vZEnd = find(diff(vZero)>1);
        vZBeg = vZEnd+1;
        vZEnd = vertcat(vZEnd,length(vZero));
        vZBeg = vertcat(1,vZBeg);

        vZShort = [];
        for i=1:length(vZBeg)
            if vZEnd(i) - vZBeg(i) < tStop
                vZShort = horzcat(vZShort,[vZBeg(i):vZEnd(i)]);
            end
        end   
        vZero(vZShort) = [];
        
        figure;
        scatter(ROISig(1+(1:length(vRot))),vRotConv(1:length(vRot)));
        
        CCNow = corrcoef(ROISig(1+(1:length(vRot))),vRotConv(1:length(vRot)));
%         CCNow = corrcoef(ROISig(1+setdiff(1:length(vRot),vZero)),vRotConv(setdiff(1:length(vRot),vZero)));
        CCs{flyID}.dark(darkID) = CCNow(2,1);
    end
    
end

NOCCs = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,3,1);
hold on;
for flyID = 1:length(allFlyData)
    scatColor = ones(length(CCs{flyID}.dark),1)*[0 0 0];
    scatter(zeros(length(CCs{flyID}.dark),1)+flyID,CCs{flyID}.dark,25,scatColor,'filled');
    line([flyID-0.4 flyID+0.4], [median(CCs{flyID}.dark) median(CCs{flyID}.dark)], 'Color','k','LineWidth',3); 
end
set(gca,'FontSize',16);
xlabel('Fly #');
ylabel('correlation (\DeltaF/F vs. v_{R})');
xlim([ 0 (numFlies+1)]);
ylim([-1 1]);
title('dark');
line([0 (numFlies+1)], [0 0],'Color','k');