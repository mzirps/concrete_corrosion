function Vf = shiftRotateCyl(Points_corr_surf_temp,i,deltax,deltay,deltaz,rebarCenter,cylHeight)

% MZ shift rebar
Points_corr_surf_shift = [Points_corr_surf_temp(:,1) + rebarCenter(i,1)...
    Points_corr_surf_temp(:,2) + rebarCenter(i,2) Points_corr_surf_temp(:,3) + rebarCenter(i,3) + cylHeight/2];

% MZ rotate rebar
theta_xz = asin(deltax(i,1)/deltaz(i,1));
theta_yz = asin(deltay(i,1)/deltaz(i,1));
Vf = [Points_corr_surf_shift(:,1) + Points_corr_surf_shift(:,3)*sin(theta_xz)...
    Points_corr_surf_shift(:,2) + Points_corr_surf_shift(:,3)*sin(theta_yz)...
    -Points_corr_surf_shift(:,3) + Points_corr_surf_shift(:,3)*cos(theta_xz) + Points_corr_surf_shift(:,3)*cos(theta_yz)];