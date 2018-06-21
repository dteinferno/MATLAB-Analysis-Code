%Calculate forward and rotational positions from dx, dy values

function [posRot, posFor, posLat] = PositionConverter(positionDat)
    Cam1RotCalibfact = 1.15;
    Cam2RotCalibfact = 0.93;
    Cam1PosCalibfact = 125;
    Cam2PosCalibfact = 153;
    flyAng = 30*pi/180;

    dxmod0 = double(positionDat.dx0).*cos(flyAng) + double(positionDat.dy0).*sin(flyAng);
    dymod0 = double(positionDat.dy0).*cos(flyAng) + double(positionDat.dx0).*sin(flyAng);
    dxmod1 = double(positionDat.dx1).*cos(flyAng) + double(positionDat.dy1).*sin(-flyAng);
    dymod1 = double(positionDat.dy1).*cos(flyAng) + double(positionDat.dx1).*sin(-flyAng);

	deltaFor = (dymod0 / Cam1PosCalibfact + dymod1 / Cam2PosCalibfact)*sqrt(2) / 2;
    deltaSide = (dymod0 / Cam1PosCalibfact - dymod1 / Cam2PosCalibfact)*sqrt(2) / 2;
    posRot(1) = positionDat.OffsetRot(1);
    RotOffset = positionDat.OffsetRot(1) - (dxmod0(1)*Cam1RotCalibfact + dxmod1(1)*Cam2RotCalibfact)/2;
    posFor(1) = positionDat.OffsetFor(1);
    posLat(1) = positionDat.OffsetLat(1);
    for i = 2:length(dxmod0)
        posRot(i) = (dxmod0(i)*Cam1RotCalibfact + dxmod1(i)*Cam2RotCalibfact)/2 + RotOffset; 
        posFor(i) = (deltaFor(i)-deltaFor(i-1))*cos(posRot(i)*pi/180) + (deltaSide(i)-deltaSide(i-1))*sin(posRot(i)*pi/180) + posFor(i-1);
        posLat(i) = (deltaFor(i)-deltaFor(i-1))*sin(posRot(i)*pi/180) - (deltaSide(i)-deltaSide(i-1))*cos(posRot(i)*pi/180) + posLat(i-1);
    end
    posRot = mod(posRot+180, 360)-180;
    
%     figure;
%     subplot(3,1,1);
%     plot(positionDat.t,positionDat.OffsetRot,'k');
%     hold on;
%     plot(positionDat.t,posRot,'r');
%     subplot(3,1,2);
%     plot(positionDat.t,positionDat.OffsetFor,'k');
%     hold on;
%     plot(positionDat.t,posFor,'r');
%     subplot(3,1,3);
%     plot(positionDat.t,positionDat.OffsetLat,'k');
%     hold on;
%     plot(positionDat.t,posLat,'r');
end