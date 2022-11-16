function [point_steel_surf,corr2steel] = moveCorrPoint(Points_corr_surf,Points_steel_surf,...
            ind_rad_steel,rebarCenterVec,non_uni_corr_layer_comb)
        
        % MZ sort points so that they match in steel a corr
        [~,ind_corr] = sort(mean(Points_corr_surf,2));
        [~,ind_steel] = sort(mean(Points_steel_surf(ind_rad_steel,:),2));
        corr2steel(ind_corr) = ind_steel;
        steel2corr(ind_steel) = ind_corr;
        
        % MZ alter steel points around radius to match corrosion pattern
        Con_quad_points(:,1:3) = Points_steel_surf;
        [theta_steel, rad_steel] = cart2pol(Con_quad_points(ind_rad_steel,1)-rebarCenterVec(corr2steel,1),...
            Con_quad_points(ind_rad_steel,2)-rebarCenterVec(corr2steel,2));
        rad_steel_adjust = rad_steel - non_uni_corr_layer_comb;
        Con_quad_points(:,4:5) = Con_quad_points(:,1:2);
        Con_quad_points(ind_rad_steel,4) = rad_steel_adjust.*cos(theta_steel)+rebarCenterVec(corr2steel,1);
        Con_quad_points(ind_rad_steel,5) = rad_steel_adjust.*sin(theta_steel)+rebarCenterVec(corr2steel,2);
        
        point_steel_surf = [Con_quad_points(:,4) Con_quad_points(:,5) Con_quad_points(:,3)];