%% Load all of the data
% Specify directories
dirRoot = 'D:\Imaging\2Color';

imDir{1} = '20160423-Flip';
imDir{2} = '20160425-Flip';
imDir{3} = '20160426-Flip';
imDir{4} = '20160429';
imDir{5} = '20160505';
redLine = '60D05';
greenLine = '37F06';

cd(strcat(dirRoot,'\',imDir{1}));

% Get the fly IDs
allFlyData{1}.ID = 'Empty';
numFlies = 1;
for dirs = 1:length(imDir)
    fileNames = dir(strcat(dirRoot,'\',imDir{dirs}));
    for fileID = 3:length(fileNames)
        fileName = fileNames(fileID).name;
        if strcmpi(fileName(end-8:end),'NORaw.mat')
            flyName = strcat(imDir{dirs},'-',fileName(1:4));
            for flyStep = 1:length(allFlyData)
                if strcmp(flyName,allFlyData{flyStep}.ID)
                    break;
                end
                if flyStep == length(allFlyData)
                    allFlyData{numFlies}.ID = flyName;
                    numFlies = numFlies + 1;
                end
            end
        end
    end
end
numFlies = numFlies-1;

% Get the position data, ROI data, and trial ID for all of the flies and all of the trials
for dirs = 1:length(imDir)
    fileNames = dir(strcat(dirRoot,'\',imDir{dirs}));
    for fileID = 3:length(fileNames)
        fileName = fileNames(fileID).name;
        if strcmpi(fileName(end-8:end),'NORaw.mat')
            flyName = strcat(imDir{dirs},'-',fileName(1:4));
            load(strcat(dirRoot,'\',imDir{dirs},'\',fileName));
            
            for flyID = 1:length(allFlyData)
                if strcmp(flyName,allFlyData{flyID}.ID)
                    fileParts = strsplit(fileName,'_');
                     if strcmpi(fileParts(5),'1x')
                        allFlyData{flyID}.gainOne{str2num(fileName(end-10))}.RROIave = RROIaveREF;
                        allFlyData{flyID}.gainOne{str2num(fileName(end-10))}.GROIave = GROIaveREF;
                    elseif strcmpi(fileParts(5),'2x')
                        allFlyData{flyID}.gainTwo{str2num(fileName(end-10))}.RROIave = RROIaveREF;
                        allFlyData{flyID}.gainTwo{str2num(fileName(end-10))}.GROIave = GROIaveREF;
                    elseif strcmpi(fileParts(5),'Dark')
                        allFlyData{flyID}.dark{str2num(fileName(end-10))}.RROIave = RROIaveREF;
                        allFlyData{flyID}.dark{str2num(fileName(end-10))}.GROIave = GROIaveREF;
                    end
                end
            end
        end 
    end
end

%% Plot R vs G
RGComp = figure('units','normalized','outerposition',[0 0 0.5 1]);

for flyID = 1:length(allFlyData)
    % Plot parameters for the dark data
    for darkID = 1:length(allFlyData{flyID}.dark);
        subplot(4,3,flyID);
        hold on;
        scatter(allFlyData{flyID}.dark{darkID}.GROIave(1,:),allFlyData{flyID}.dark{darkID}.RROIave(1,:),10,'filled');
        scatter(allFlyData{flyID}.dark{darkID}.GROIave(2,:),allFlyData{flyID}.dark{darkID}.RROIave(2,:),10,'filled');
        set(gca,'FontSize',14);
        xlabel('Fgreen');
        ylabel('Fred');
        axis tight;
        title(strcat('Fly #',num2str(flyID)));
    end
end

set(RGComp,'PaperPositionMode','manual','PaperOrientation','portrait','PaperUnits','inches','PaperPosition',[0 0 8.5 11]);
print(RGComp,'..\Results\RvsGComp','-dpdf');