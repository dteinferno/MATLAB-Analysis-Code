% 11/22/2014
% Imaging Data Analysis Script for 60D05/translation experiments
% Dan Turner-Evans

%% Plot summary data from all trials
% Load the position information
allPathname = uigetdir('D:/Imaging','Select the directory');
fileNames = dir(allPathname);
flyNum = input('Fly #?');
allPaths = figure('units','normalized','outerposition',[0 0 1 1])
spInc = 1;
for fID = 3:length(fileNames)
    fName = fileNames(fID).name; 
    if (strcmp(fName(1:3),'Fly') & strcmp(fName(4),num2str(flyNum)) & ~strcmp(fName(end-7:end-4),'SYNC') & strcmpi(fName(end-2:end),'txt'))
        fileID = fopen(strcat(allPathname,'\',fName));
        tstamp = fgetl(fileID);
        exTypeStr = strsplit(tstamp,'_');
        exType = exTypeStr{end}(1:end-4);
        formatSpec = '%s %f %s %f %s %f %s %f %s %f %s %f %s %d %s %d %s %d %s %d %s %d %s %f';
        N=400000;
        C = textscan(fileID,formatSpec,N,'CommentStyle','Current','Delimiter','\t');
        t = C{1,2}; % Time
        OffsetRot = C{1,4}; % Stripe rotational offset
        OffsetRot = mod(OffsetRot+180, 360)-180;
        OffsetFor = C{1,6}; % Stripe forward offset
        OffsetLat = C{1,8}; % Stripe lateral offset
        dx0 = C{1,10}; % X position of the ball from camera 1 
        dx1 = C{1,12}; % X position of the ball from camera 2
        dy0 = C{1,14};
        dy1 = C{1,16};
        closed = C{1,18};
        direction = C{1,20};
        trans = C{1,22};
        gain = C{1,24};
        fclose(fileID);

        % Discriminate between the different trials
        trialbreaks1 = find(abs(diff(closed))>=1);
        trialbreaks2 = find(abs(diff(direction))>=1);
        trialbreaks = union(trialbreaks1,trialbreaks2);
        trialbreaks = vertcat(0,trialbreaks,length(t));
        numTrials = length(trialbreaks)-1;

        % Downsample the data to 20Hz
        sampFreq = 1/mean(diff(t));
        pts2sum = round(sampFreq/20);
        tDS = zeros(floor(length(OffsetLat)/pts2sum),1);
        OffsetRotDS = tDS;
        OffsetForDS = tDS;
        OffsetLatDS = tDS;
        for DSstep = 1:length(tDS);
            tDS(DSstep) = mean(t(1+pts2sum*(DSstep-1):pts2sum*DSstep));
            OffsetRotDS(DSstep) = mean(OffsetRot(1+pts2sum*(DSstep-1):pts2sum*DSstep));
            OffsetForDS(DSstep) = mean(OffsetFor(1+pts2sum*(DSstep-1):pts2sum*DSstep));
            OffsetLatDS(DSstep) = mean(OffsetLat(1+pts2sum*(DSstep-1):pts2sum*DSstep));
        end
        
        % Plot the fly's position in the arena throughout the trial
        subplot(3,4,spInc);
        hold on;
        c = colormap(hsv(length(tDS)));
        scatter(-OffsetLatDS,OffsetForDS,1,c);
        if (strcmp(exType,'LightCylFlatBack') || strcmp(exType,'LightCyl3ObjBack') || strcmp(exType,'LightCylFlatBackInv'))
            viscircles([0 0], 0.5,'EdgeColor','b');
        end
        if strcmp(exType,'LightCylAlone') || strcmp(exType,'LightCylAloneInv')
            viscircles([0 0], 0.5,'EdgeColor','b');
            rectangle('Position',[-0.5 -0.5 1 1],'Curvature',1,'FaceColor','k');
        end
        if (strcmp(exType,'NoCyl3ObjBack') || strcmp(exType,'LightCyl3ObjBack'))
            line([0 2*cos(pi/6)],[0 2*sin(pi/6)],'LineStyle','--','Color',[0 0 0]);
            line([0 -2*cos(pi/6)],[0 2*sin(pi/6)],'LineStyle','--','Color',[0 1 0]);
            line([0 2*cos(pi/6-pi/2)],[0 2*sin(pi/6-pi/2)],'LineStyle','--','Color',[1 0 0]);
        end
        axis equal;
        
        fParts = strsplit(fName,'_');
        title(strcat(exType,'-',fName(end-5:end-4)),'FontSize',10);
        xlabel('Distance (cm)');
        ylabel('Distance (cm)');
        set(gca,'FontSize',10);
        scatter(0,-7,100,'p','k')
        
        spInc = spInc+1;
        
        % Plot the statistics for each fly
        flyStats = figure('units','normalized','outerposition',[0 0 1 1]);

        % Show the fly's path
        subplot(2,2,1);
        hold on;
        scatter(-OffsetLatDS,OffsetForDS,1,c);
        for plotPt = 1:length(tDS)
           line([-OffsetLatDS(plotPt) -OffsetLatDS(plotPt) - 0.15*sin(OffsetRotDS(plotPt)*pi/180)],[OffsetForDS(plotPt) OffsetForDS(plotPt) + 0.15*cos(OffsetRotDS(plotPt)*pi/180)], 'Color', c(plotPt,:)); 
        end
        if (strcmp(exType,'LightCylFlatBack') || strcmp(exType,'LightCyl3ObjBack') || strcmp(exType,'LightCylFlatBackInv'))
            viscircles([0 0], 0.5,'EdgeColor','b');
        end
        if strcmp(exType,'LightCylAlone') || strcmp(exType,'LightCylAloneInv')
            viscircles([0 0], 0.5,'EdgeColor','b');
            rectangle('Position',[-0.5 -0.5 1 1],'Curvature',1,'FaceColor','k');
        end
        if (strcmp(exType,'NoCyl3ObjBack') || strcmp(exType,'LightCyl3ObjBack'))
            line([0 2*cos(pi/6)],[0 2*sin(pi/6)],'LineStyle','--','Color',[0 0 0]);
            line([0 -2*cos(pi/6)],[0 2*sin(pi/6)],'LineStyle','--','Color',[0 1 0]);
            line([0 2*cos(pi/6-pi/2)],[0 2*sin(pi/6-pi/2)],'LineStyle','--','Color',[1 0 0]);
        end
        axis equal;
        scatter(0,-7,100,'p','k')
        fParts = strsplit(fName,'_');
        title(strcat(exType,'-',fName(end-5:end-4)),'FontSize',16);
        xlabel('Distance (cm)');
        ylabel('Distance (cm)');
        set(gca,'FontSize',16);

        % Show the fly's forward velocity distribution
        subplot(2,2,3);
        frameRate = 1/mean(diff(tDS));
        netSpeed = sqrt(diff(OffsetForDS).^2+diff(OffsetLatDS).^2)*frameRate;
        [vCounts, vVals] = hist(netSpeed,round(frameRate*10));
        scatter(vVals,vCounts,'filled');
        hold on;
%         pks = findpeaks(smooth(vCounts,5));
%         pkMax = find(smooth(vCounts,5)==max(pks));
%         xOver = min(find(vCounts==min(vCounts(1:pkMax))));
%         if xOver>3
            ft1 = fittype('exp1');
%             cf1 = fit(vVals(2:xOver-1)',vCounts(2:xOver-1)',ft1);
%             plot(cf1,vVals(2:xOver-1)',vCounts(2:xOver-1)');
%         end
        ft2 = fittype('gauss1');
%         cf2 = fit(vVals(xOver+2:end)',vCounts(xOver+2:end)',ft2);
%         plot(cf2,vVals(xOver+2:end)',vCounts(xOver+2:end)');
        xlabel('Forward speed (cm/sec)');
        ylabel('Counts');
        axis tight;
        ylim([0 max(vCounts(2:end))]);
        xlim([0 2.5]);
        set(gca,'FontSize',16);
%         line([vVals(xOver) vVals(xOver)], [0 max(vCounts(2:end))],'LineStyle','--','Color',[0 0 0],'LineWidth',2);
%         legend off;
%         perSlow = 100*sum(vCounts(1:xOver))/sum(vCounts);
%         text(vVals(round(xOver/2)),max(vCounts(2:end))/2,strcat(num2str(round(perSlow)),'%'),'FontSize',14);
%         perFast = 100*sum(vCounts(xOver+1:end))/sum(vCounts);
%         text(vVals(min(pkMax)),max(vCounts(2:end))/2,strcat(num2str(round(perFast)),'%'),'FontSize',14);

        % Show the fly's rotational velocity distribution
        subplot(2,2,4);
        rotSpeed = diff(OffsetRotDS)*frameRate*pi/180;
        [rCounts, rVals] = hist(rotSpeed,round(frameRate*60));
        scatter(rVals,rCounts,'filled');
        hold on;
        cf3 = fit(rVals(find(rVals>pi/16 & rVals < 2*pi))',rCounts(find(rVals>pi/16 & rVals < 2*pi))',ft1);
        cf4 = fit(rVals(find(rVals<-pi/16 & rVals > -2*pi))',rCounts(find(rVals<-pi/16 & rVals > -2*pi))',ft1);
        plot(cf3,rVals(find(rVals>pi/16 & rVals < 2*pi))',rCounts(find(rVals>pi/16 & rVals < 2*pi))');
        plot(cf4,rVals(find(rVals<-pi/16 & rVals > -2*pi))',rCounts(find(rVals<-pi/16 & rVals > -2*pi))');
        xlabel('Rotational velocity (rad/sec)');
        ylabel('Counts');
        axis tight;
        xlim([-pi pi]);
        set(gca,'FontSize',16);
        legend off;
        line([0 0], [0 max(rCounts)],'LineStyle','--','Color',[0 0 0],'LineWidth',2);
        perRight = -100*(cf3.a/cf3.b)/(-cf3.a/cf3.b+cf4.a/cf4.b);
        text(pi/4,max(rCounts)/2,strcat(num2str(round(perRight)),'%'),'FontSize',14);
        perLeft = 100*(cf4.a/cf4.b)/(-cf3.a/cf3.b+cf4.a/cf4.b);
        text(-pi/4,max(rCounts)/2,strcat(num2str(round(perLeft)),'%'),'FontSize',14);

        % Show the discrepancies between the fly's body heading and rotational
        % heading
%         subplot(3,6,[16:18]);
%         latDiff = diff(OffsetLatDS);
%         forDiff = diff(OffsetForDS);
%         bodHead = atan2(latDiff,forDiff);
%         for i=2:length(bodHead)
%             if latDiff(i) < 2E-4 & forDiff < 2E-4
%                 bodHead(i) = bodHead(i-1);
%             end
%         end
%         OffsetRotUnwrap = OffsetRotDS*pi/180;
%         [headCounts, headVals] = hist(bodHead-OffsetRotUnwrap(2:length(bodHead)+1),round(frameRate*10));
%         scatter(headVals,headCounts,'filled');
%         hold on;
%         cf5 = fit(headVals(find(abs(headVals) < pi))',headCounts(find(abs(headVals) < pi))',ft2);
%         plot(cf5,headVals(find(abs(headVals) < pi))',headCounts(find(abs(headVals) < pi))');
%         xlabel('Body angle - heading deviation (rad per 1 sec)');
%         ylabel('Counts');
%         axis tight;
%         xlim([-2*pi 2*pi]);
%         set(gca,'FontSize',16);
%         legend off;
%         line([cf5.b1 cf5.b1], [0 max(headCounts)],'LineStyle','--','Color',[0 0 0],'LineWidth',2);
%         fit95 = confint(cf5);
%         text(-pi,max(headCounts)/2,strcat('A: ',num2str(cf5.a1),'(+/-',num2str(round(100*(cf5.a1-fit95(1))/cf5.a1)),'%)'),'FontSize',14);
%         text(-pi,max(headCounts)/3,strcat('x_{0}: ',num2str(cf5.b1),'(+/-',num2str(round(100*(cf5.b1-fit95(3))/cf5.b1)),'%)'),'FontSize',14);
%         text(-pi,max(headCounts)/6,strcat('\sigma: ',num2str(cf5.c1),'(+/-',num2str(round(100*(cf5.c1-fit95(5))/cf5.c1)),'%)'),'FontSize',14);
        print(flyStats,strcat(allPathname,'\',fParts{1},'_Trial_',fName(end-5:end-4),'_Stats'),'-djpeg');
        set(flyStats,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
        print(flyStats,strcat(allPathname,'\',fParts{1},'_Trial_',fName(end-5:end-4),'_Stats'),'-dpdf');
        delete(flyStats);
    end
end
print(allPaths,strcat(allPathname,'\',fParts{1},'_Traj'),'-djpeg');
set(allPaths,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
print(allPaths,strcat(allPathname,'\',fParts{1},'_Traj'),'-dpdf');
delete(allPaths);
% 
% 
% % Write html code to make it easier to look at later
% fileID = fopen(strcat(allPathname,'\',fParts{1},'.htm'),'w');
% fprintf(fileID,'<html>\n');
% fprintf(fileID,'<head>\n<meta http-equiv="Content-Type" content="text/html; charset=windows-1252">\n');
% fprintf(fileID,'<title>%s</title>\n',fParts{1});
% fprintf(fileID,'<head>\n');
% fprintf(fileID,'<body topmargin="0" leftmargin="0" rightmargin="0" bottommargin="0" marginheight="0" marginwidth="0" bgcolor="#FFFFFF">\n');
% for fID = 3:length(fileNames)
%     fName = fileNames(fID).name; 
%     if (strcmp(fName(1:3),'Fly') & strcmp(fName(4),num2str(flyNum)) & ~strcmp(fName(end-7:end-4),'SYNC') & strcmpi(fName(end-2:end),'txt'))
%         fParts = strsplit(fName,'_');
%         fprintf(fileID,'<A HREF="%s">\n',strcat(fName(1:end-4),'.htm'));
%         fprintf(fileID,'<IMG SRC="%s" ALT="%s" WIDTH=550 HEIGHT=425>\n',strcat(fParts{1},'_Trial',fName(end-5:end-4),'_Stats.jpg'),strcat(fParts{1},'_Trial',fName(end-5:end-4)));
%         fprintf(fileID,'</A>\n');
%     end
% end
% fprintf(fileID,'</body>\n');
% fprintf(fileID,'</html>\n');
% fclose(fileID);