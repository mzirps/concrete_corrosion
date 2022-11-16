function [Corr_layer_thick_max,F_diff_non_tot,rebarCenterVec,non_uni_corr_layer] = findCorrDim(num_longRebar,...
            Points_corr_surf,exp_z_corr_layer,exp_theta_corr_layer,rebar_point_count,...
            rebarCenter,delta_T,l_A,rebarR,kk,V_cp_min,V_cp_max,CAR_0,CAR_max,...
            coeff_por,n_non,exp_corr_layer,non_uni_corr_layer)
        
        Corr_layer_thick_max=cell(num_longRebar,1);
        F_diff_non=cell(num_longRebar,1); % force caused by corrosion
        count_rebar_surf_points = 1;
        rebarCenterVec = ones(size(Points_corr_surf,1),2);
        
        for i_rebar=1:num_longRebar
%             disp('Calculating CAR for rebar ' + string(i_rebar));
            [exp_z, exp_theta] = ndgrid(exp_z_corr_layer{i_rebar},exp_theta_corr_layer{i_rebar});
            F_non_uni_corr_layer = scatteredInterpolant(exp_theta(:),exp_z(:),exp_corr_layer{i_rebar,kk}(:),'linear','nearest');
            
            Points_i_rebar_surf = Points_corr_surf(count_rebar_surf_points:(rebar_point_count{i_rebar,1}+(count_rebar_surf_points-1)),:);
            [theta_con,~] = cart2pol(Points_i_rebar_surf(:,1)-rebarCenter(i_rebar,1),Points_i_rebar_surf(:,2)-rebarCenter(i_rebar,2));
            theta_con = rad2deg(theta_con);
            theta_con(theta_con(:,1)<0,1) = theta_con(theta_con(:,1)<0,1)+360;
            non_uni_corr_layer{i_rebar,kk} = F_non_uni_corr_layer(theta_con,Points_i_rebar_surf(:,3));
            Corr_layer_thick_max{i_rebar} = max(non_uni_corr_layer{i_rebar,kk});
            
            % MZ determine rebar radius as vector for later
            rebarCenterVec(count_rebar_surf_points:(rebar_point_count{i_rebar,1}...
                +(count_rebar_surf_points-1)),:) = [rebarCenter(i_rebar,1)...
                *rebarCenterVec(count_rebar_surf_points:(rebar_point_count{i_rebar,1}+...
                (count_rebar_surf_points-1)),1) rebarCenter(i_rebar,2)...
                *rebarCenterVec(count_rebar_surf_points:(rebar_point_count{i_rebar,1}+(count_rebar_surf_points-1)),2)];
            count_rebar_surf_points = count_rebar_surf_points + rebar_point_count{i_rebar,1};
            
            %wu non_uni_corr_layer is the corrosion layer thicknesss around
            %steel surface;
            
            
            % see ref 1 for the equation
            %eq 23
            %V_cp_non:expanded volume of corrosion products
            V_cp_non = (3.*delta_T.*pi.*l_A(i_rebar,1).*((rebarR(i_rebar)).^2-...
                (rebarR(i_rebar)-non_uni_corr_layer{i_rebar,kk}).^2))*(1000^3);
            kappa = zeros(size(V_cp_non,1),1);
            %eq 28
            kappa(V_cp_non <= V_cp_min(i_rebar,1)) = 0;
            kappa(V_cp_non > V_cp_max(i_rebar,1)) = 1;
            kappa(V_cp_min(i_rebar,1) < V_cp_non & V_cp_non <= V_cp_max(i_rebar,1)) = ...
                (V_cp_non(V_cp_min(i_rebar,1) < V_cp_non & V_cp_non <= V_cp_max(i_rebar,1))...
                - V_cp_min(i_rebar,1))/(V_cp_max(i_rebar,1)-V_cp_min(i_rebar,1));
            
            %eq 27
            CAR = CAR_0+(CAR_max-CAR_0).*kappa;
            
            %eq26
            V_cm = pi.*l_A(i_rebar,1)*1000.*((rebarR(i_rebar)*1000+CAR).^2-(rebarR(i_rebar)*1000).^2);
            
            %eq 24
            V_CAR = coeff_por.*V_cm; %wu coeff_por: capillary porosity
            
            % lambda_CAR_non: penetration of corrosion products into the accessible pore space of the cementitious material
            lambda_CAR_non = zeros(size(V_cp_non,1),1);
            %eq 22
            lambda_CAR_non(V_cp_non <  V_CAR) = (V_cp_non(V_cp_non <  V_CAR)./V_CAR(V_cp_non <  V_CAR)).^n_non;
            lambda_CAR_non(V_cp_non >= V_CAR) = 1;
            
            F_diff_non{i_rebar,1} = lambda_CAR_non.*delta_T; %wu adjusted temperature increment eq.21
            
        end
        
        F_diff_non_tot = cat(1, F_diff_non{:});