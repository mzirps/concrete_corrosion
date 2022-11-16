function [Points_corr_surf,Points_steel_surf,F_corr_surf,F_steel_surf,ind_rad_steel,...
            rebar_point_count,steel_point_count,rebar_F_count,steel_F_count,...
            top,bottom] = formSteelCorrMesh(num_longRebar,rebarR,No_ele_rebar_surf,...
            No_ele_rebar_vert,deltax,deltay,deltaz,rebarCenter)
        
        Points_corr_surf = [];
        Points_steel_surf = [];
        F_corr_surf = [];
        F_steel_surf = [];
        ind_rad_steel = [];
        count_ind = 0;
        rebar_point_count = cell(num_longRebar,1);
        steel_point_count = cell(num_longRebar,1);
        rebar_F_count = cell(num_longRebar,1);
        steel_F_count = cell(num_longRebar,1);
        top = cell(num_longRebar,1);
        bottom = cell(num_longRebar,1);
        for i = 1:num_longRebar
%             disp('patching rebar ' + string(i) + ' corrosion surface')
            % MZ define dimensions
            input_corr.cylRadius = rebarR(i,1);
            input_corr.numRadial = No_ele_rebar_surf;
            input_corr.numHeight = No_ele_rebar_vert;
            input_corr.cylHeight = sqrt(deltax(i,1)^2 + deltay(i,1)^2 + deltaz(i,1)^2);
            input_corr.meshType = 'tri';
            input_corr.closeOpt = 0;
            
            % MZ input into function
            [F_corr_surf_temp, Points_corr_surf_temp] = patchcylinder(input_corr);
            
            % MZ alter inputs to steel surf
            input_corr.closeOpt = 1;
            
            % MZ create steel surf
            [F_steel_surf_temp, Points_steel_surf_temp] = patchcylinder(input_corr);
            
            % MZ get outer indexing of steel
            [~,rho] = cart2pol(Points_steel_surf_temp(:,1),Points_steel_surf_temp(:,2));
            ind_rad_steel = [ind_rad_steel; count_ind + find(rho > rebarR(i,1)-.0001)];
            count_ind = count_ind + size(Points_steel_surf_temp,1);
            
            % MZ need to shift and rotate rebar linearly to get into
            % correct location
            Points_corr_surf_shift_rot = shiftRotateCyl(Points_corr_surf_temp,i,deltax,deltay,deltaz,rebarCenter,input_corr.cylHeight);
            Points_steel_surf_shift_rot = shiftRotateCyl(Points_steel_surf_temp,i,deltax,deltay,deltaz,rebarCenter,input_corr.cylHeight);
            
            
            % MZ count of points in each reinforcement
            rebar_point_count{i,1} = size(Points_corr_surf_shift_rot,1);
            steel_point_count{i,1} = size(Points_steel_surf_shift_rot,1);
            rebar_F_count{i,1} = size(F_corr_surf_temp,1);
            steel_F_count{i,1} = size(F_steel_surf_temp,1);
            
            % MZ define top and bottom of rebar
            top{i,1} = max(Points_steel_surf_shift_rot(:,3));
            bottom{i,1} = min(Points_steel_surf_shift_rot(:,3));
            
            % MZ add rebar to overall matrix
            F_corr_surf = [F_corr_surf; F_corr_surf_temp + size(Points_corr_surf,1)];
            F_steel_surf = [F_steel_surf; F_steel_surf_temp + size(Points_steel_surf,1)];
            Points_corr_surf = [Points_corr_surf; Points_corr_surf_shift_rot];
            Points_steel_surf = [Points_steel_surf; Points_steel_surf_shift_rot];
        end