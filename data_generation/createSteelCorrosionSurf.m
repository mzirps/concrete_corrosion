function [Vsc,Fsc,Fsrf,Fcf,Ics] = createSteelCorrosionSurf(Vs,Fs,Vc,Fc,t,b,i_rad_steel,centS,centC)

% find indexing of top and bottom radius (keep separate)
indexTopS = i_rad_steel(Vs(i_rad_steel,3) == t);
indexBotS = i_rad_steel(Vs(i_rad_steel,3) == b);

indexTopC = find(Vc(:,3) == t);
indexBotC = find(Vc(:,3) == b);

% move points to polar coordinates
[thetaTopS, ~] = cart2pol(Vs(indexTopS,1)-centS(Vs(i_rad_steel,3) == t,1),Vs(indexTopS,2)-centS(Vs(i_rad_steel,3) == t,2));
[thetaBotS, ~] = cart2pol(Vs(indexBotS,1)-centS(Vs(i_rad_steel,3) == b,1),Vs(indexBotS,2)-centS(Vs(i_rad_steel,3) == b,2));

[thetaTopC, ~] = cart2pol(Vc(indexTopC,1)-centC(indexTopC,1),Vc(indexTopC,2)-centC(indexTopC,2));
[thetaBotC, ~] = cart2pol(Vc(indexBotC,1)-centC(indexBotC,1),Vc(indexBotC,2)-centC(indexBotC,2));

indexTopS = indexTopS + size(Vc,1);
indexBotS = indexBotS + size(Vc,1);

% sort points so that indexing is in order of circumference
[~,indexTopSSort] = sort(thetaTopS);
[~,indexBotSSort] = sort(thetaBotS);

[~,indexTopCSort] = sort(thetaTopC);
[~,indexBotCSort] = sort(thetaBotC);

% rearrange indexing
indexTopSSortShift = indexTopS(indexTopSSort,1);
indexBotSSortShift = indexBotS(indexBotSSort,1);

indexTopCSortShift = indexTopC(indexTopCSort,1);
indexBotCSortShift = indexBotC(indexBotCSort,1);

% create indexing for cap between cylinders
indexTopSShift = [indexTopSSortShift(2:end,1); indexTopSSortShift(1,1)];
indexTopCShift = [indexTopCSortShift(end,1); indexTopCSortShift(1:(end-1),1)];
    
indexBotSShift = [indexBotSSortShift(2:end,1); indexBotSSortShift(1,1)];
indexBotCShift = [indexBotCSortShift(end,1); indexBotCSortShift(1:(end-1),1)];

Fcap = [indexTopCSortShift indexTopSShift indexTopSSortShift ;
    indexTopCSortShift indexTopCShift indexTopSSortShift;
    indexBotCSortShift indexBotSShift indexBotSSortShift;
    indexBotCSortShift indexBotCShift indexBotSSortShift];

% find indexign of caps of steel
indicator = zeros(size(Vs,1),1);
indicator(Vs(:,3) > t-.000001,1) = 1;
indicator(Vs(:,3) < b+.000001,1) = 1;

FsCap = Fs(mean(indicator(Fs),2) == 1,:);

FsCap = FsCap + size(Vc,1);

% find Fsrf
Fsrf = findFsrf(Vs,Fs,Vc,Fc,i_rad_steel,centS,centC);
  
% combine points to create Vsc and Fsc
Vsc = [Vc; Vs];
Fsc = [Fc; Fcap; FsCap];

% create indexing of corrosion surface
Fcf = [Fcap;Fc;Fsrf];

% Find indexing of corrosion surface
Ics = [(1:size(Vc,1))'; i_rad_steel+size(Vc,1)];
save('index_check_2')
