%% Load the imaging and DAQ data from a 2D VR experiment
function [fullpath, RstackMaxInt, GstackMaxInt, RstackMean, GstackMean, positionDat] = VRImDatLoadDiscardNoPrompt2Color(imageFilename,imagePathname,posFilename,posPathname,numFiles,planesDiscard)

% Get the 2P info;
fullpath = strcat(imagePathname,imageFilename);
cd(imagePathname);

fullpath = {};
num_images = 0;
num_images_ind = zeros(numFiles,1);
info = {};
for fID = 1:numFiles
    fullpath{fID} = strcat(imagePathname,imageFilename(1:end-4),'.tif');
    info{fID} = imfinfo(fullpath{fID});
    num_images_ind(fID) = numel(info{fID});
    num_images = num_images + num_images_ind(fID);
end
recSpecs = info{1}(1).ImageDescription;
planeLoc = strfind(recSpecs, 'numFramesPerVolume');
discardLoc = strfind(recSpecs, 'numDiscardFlybackFrames');
num_planes = str2num(recSpecs(planeLoc+21:planeLoc+22));
num_discards = str2num(recSpecs(discardLoc+26));
if isempty(num_planes)
    hTiff = Tiff(strcat(imagePathname,imageFilename));
    configStr = hTiff.getTag('Software');
    planeLoc = strfind(configStr, 'numFramesPerVolume');
    discardLoc = strfind(configStr, 'numDiscardFlybackFrames');
    num_planes = str2num(configStr(planeLoc+21:planeLoc+22));
    num_discards = str2num(configStr(discardLoc+26));
end

% Read in the stack
% Get the file info and # of planes and frames
width = info{1}(1).Width;
height = info{1}(1).Height;
numFrames = floor(num_images/num_planes);
num_images_ind(end) = num_images_ind(end) - (num_images - num_planes*numFrames); 

% Find the mean volume projection
RstackMean = double(zeros(height,width,numFrames/2));
GstackMean = double(zeros(height,width,numFrames/2));
h = waitbar(0.0,'Loading TIFF stack (mean calc)...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Loading TIFF stack...');
imOffset = 0;
for fID=1:numFiles
    for incIm = 1:num_images_ind(fID)
        if mod(incIm+imOffset,100)==0
            waitbar((incIm+imOffset)/num_images,h,['Loading frame# ' num2str(incIm+imOffset) ' out of ' num2str(num_images)]);
        end
        if mod(ceil((incIm+imOffset)/2),num_planes) ~= 0 & mod(ceil((incIm+imOffset)/2),num_planes) ~= 0
            if mod(incIm+imOffset,2)
                RstackMean(:,:,ceil((incIm+imOffset)/(2*num_planes))) = double(imread(fullpath{fID}, incIm, 'Info', info{fID}))+ RstackMean(:,:,ceil((incIm+imOffset)/(2*num_planes)));
            else
                GstackMean(:,:,ceil((incIm+imOffset)/(2*num_planes))) = double(imread(fullpath{fID}, incIm, 'Info', info{fID}))+ GstackMean(:,:,ceil((incIm+imOffset)/(2*num_planes)));
            end
        end
    end
    imOffset = imOffset+numel(info{fID});
end
delete(h);

% Find the max intensity projection
RstackMaxInt = double(zeros(height,width,numFrames/2));
GstackMaxInt = double(zeros(height,width,numFrames/2));
RminiStack = double(zeros(height,width,(num_planes-num_discards)));
GminiStack = double(zeros(height,width,(num_planes-num_discards)));
h = waitbar(0.0,'Loading TIFF stack (mean calc)...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Loading TIFF stack...');
imOffset = 0;
for fID=1:numFiles
    for incIm = 1:num_images_ind(fID)
        if mod(incIm+imOffset,100)==0
            waitbar((incIm+imOffset)/num_images,h,['Loading frame# ' num2str(incIm+imOffset) ' out of ' num2str(num_images)]);
        end
        if mod(ceil((incIm+imOffset)/2),num_planes) ~= 0
            if mod(incIm+imOffset,2)
                RminiStack(:,:,mod(ceil((incIm+imOffset)/2),num_planes)) = double(imread(fullpath{fID}, incIm, 'Info', info{fID}));
            else
                GminiStack(:,:,mod(ceil((incIm+imOffset)/2),num_planes)) = double(imread(fullpath{fID}, incIm, 'Info', info{fID}));
            end
            if (mod((incIm+imOffset)/2,num_planes) == num_planes-num_discards)
                RstackMaxInt(:,:,ceil((incIm+imOffset)/(2*num_planes))) = max(RminiStack,[],3);
                GstackMaxInt(:,:,ceil((incIm+imOffset)/(2*num_planes))) = max(GminiStack,[],3);
            end
        end
    end
    imOffset = imOffset+numel(info{fID});
end
delete(h);

% Get the DAQ info
SYNCFilename = strcat(posFilename(1:end-4),'_SYNC',posFilename(end-3:end));
SYNCPathname = posPathname;

% Pull out the fly positional information
fileID = fopen(strcat(posPathname,posFilename));
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

numDatPts = length(OffsetFor);

positionDat = {};
positionDat.t = t(1:numDatPts);
positionDat.OffsetRot = OffsetRot(1:numDatPts);
positionDat.OffsetFor = OffsetFor;
positionDat.OffsetLat = OffsetLat;
positionDat.dx0 = dx0;
positionDat.dx1 = dx1;
positionDat.dy0 = dy0;
positionDat.dy1 = dy1;
positionDat.closed = closed;
positionDat.direction = direction;
positionDat.trans = trans;
positionDat.gain = gain;
positionDat.exType = exType;


% Load photodiode and frame time stamps
timeStamps = importdata(strcat(SYNCPathname,SYNCFilename));

% Pull out the time stamps for the frame grab signal
% clear tFrameGrab;
% sampleData = 1;
% upperLim = max(timeStamps(:,1));
% offset = round(0.8/(88)*10000);
% startFrameGrab = find(timeStamps(:,1) > 0.9*upperLim);
% incDat = startFrameGrab(1)-2;
% inct = 1;
% while (sampleData)
%     if (timeStamps(incDat+1,1) > 0.9*upperLim && timeStamps(incDat-2,1) < timeStamps(incDat+1,1))
%         tFrameGrab(inct) = incDat+1;
%         inct = inct +1;
%         incDat = incDat + offset;
%     end
%     incDat=incDat+1;
%     if incDat > length(timeStamps)-1
%         break
%     end
% end

clear tFrameGrab;
sampleData = 1;
upperLim = max(timeStamps(:,1));
startFrameGrab = find(timeStamps(:,1) > 0.8*upperLim);
incDat = startFrameGrab(1)-2;
inct = 1;
while (sampleData)
    if (timeStamps(incDat+1,1) > 0.8*upperLim && timeStamps(incDat,1) < 0.8*upperLim)
        tFrameGrab(inct) = incDat+1;
        inct = inct +1;
        incDat = incDat + 1;
    end
    incDat=incDat+1;
    if incDat > length(timeStamps)-1
        break
    end
end

positionDat.tFrameGrab = tFrameGrab;

% Pull out the time stamps for the VR refresh
clear tVR;
sampleData = 1;
upperLim = max(timeStamps(:,2));
offset = round(0.6/(360)*10000);
VRthresh = 0.8;
startVR = find(timeStamps(:,2) > VRthresh*upperLim);
incDat = startVR(1)-2;
inct = 1;
while (sampleData)
    if (timeStamps(incDat+1,2) < VRthresh*upperLim && (timeStamps(incDat-1,2) < timeStamps(incDat+1,2) || timeStamps(incDat,2) < timeStamps(incDat+1,2)))
        tVR(inct) = incDat+1;
        inct = inct +1;
        incDat = incDat + offset;
    end
    incDat=incDat+1;
    if incDat > length(timeStamps)-1
        break
    end
end

positionDat.tVR = tVR;

if 0
% Pull out the time stamps for the stimulation
clear tStim;
sampleData = 1;
upperLim = max(timeStamps(:,3));
startFrameGrab = find(timeStamps(:,3) > 0.9*upperLim);
incDat = startFrameGrab(1)-2;
inct = 1;
while (sampleData)
    if (timeStamps(incDat+1,3) > 0.9*upperLim && timeStamps(incDat,3) < 0.9*timeStamps(incDat+1,3))
        tStim(inct) = incDat+1;
        inct = inct +1;
        incDat = incDat + offset;
    end
    incDat=incDat+1;
    if incDat > length(timeStamps)-1
        break
    end
end

positionDat.tStim = tStim;

end

end
