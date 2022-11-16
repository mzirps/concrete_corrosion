function Fsrf = findFsrf(Vs,Fs,Vc,Fc,i_rad_steel,centS,centC)

% Reduce # of elements in Fs to cover just radial sections
indicator = zeros(size(Vs,1),1);
indicator(i_rad_steel,1) = 1;
Fsrf_unordered = Fs(mean(indicator(Fs),2)==1,:);

% Convert to polar coord
[thetaS,rhoS] = cart2pol(Vs(i_rad_steel,1)-centS(:,1),Vs(i_rad_steel,2)-centS(:,2));
[thetaC,rhoC] = cart2pol(Vc(:,1)-centC(:,1),Vc(:,2)-centC(:,2));

thetaS_temp = zeros(size(Vs,1),1);
rhoS_temp = zeros(size(Vs,1),1);

thetaS_temp(i_rad_steel,1) = thetaS;
rhoS_temp(i_rad_steel,1) = rhoS;

polS = [thetaS_temp rhoS_temp Vs(:,3)];
polC = [thetaC rhoC Vc(:,3)];

% Rearrange Fs to match the order of Fc
Vc_mean = mean(polC,2);
[~,ind_c] = sort(mean(Vc_mean(Fc),2));
Vs_mean = mean(polS,2);
[~,ind_s] = sort(mean(Vs_mean(Fsrf_unordered),2));
steel2corr(ind_c) = ind_s;
Fsrf = Fsrf_unordered(steel2corr,:) + size(Vc,1);