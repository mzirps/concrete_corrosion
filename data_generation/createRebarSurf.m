function [F_steelcorr_surf,F_corr_surf_final,Index_corr_surf,Index_corr_outer,...
            V_steelcorr_surf] = createRebarSurf(corr2steel,rebarCenterVec,num_longRebar,...
            point_steel_surf,steel_point_count,F_steel_surf,steel_F_count,Points_corr_surf,...
            rebar_point_count,F_corr_surf,rebar_F_count,top,bottom,ind_rad_steel)
        
        count_steel_points = 0;
        count_steel_F = 0;
        count_corr_points = 0;
        count_corr_F = 0;
        centS = rebarCenterVec(corr2steel,:);
        V_steelcorr_surf = [];
        F_steelcorr_surf = [];
        F_corr_surf_final = [];
        Index_corr_surf = [];
        Index_corr_outer = [];
        Index_corr_inner = [];
        for i = 1:num_longRebar
            [Vsc,Fsc,~,Fcf,Ics] = createSteelCorrosionSurf(point_steel_surf((count_steel_points+1):(count_steel_points+steel_point_count{i,1}),:),...
                F_steel_surf((count_steel_F+1):(count_steel_F+steel_F_count{i,1}),:)-count_steel_points,...
                Points_corr_surf((count_corr_points+1):(count_corr_points+rebar_point_count{i,1}),:),...
                F_corr_surf((count_corr_F+1):(count_corr_F+rebar_F_count{i,1}),:)-count_corr_points,...
                top{i,1},bottom{i,1},...
                ind_rad_steel((count_corr_points+1):(count_corr_points+rebar_point_count{i,1}),:)-count_steel_points,...
                centS((count_corr_points+1):(count_corr_points+rebar_point_count{i,1}),:),...
                rebarCenterVec((count_corr_points+1):(count_corr_points+rebar_point_count{i,1}),:));
            
            F_steelcorr_surf = [F_steelcorr_surf; Fsc+size(V_steelcorr_surf,1)];
            F_corr_surf_final = [F_corr_surf_final; Fcf+size(V_steelcorr_surf,1)];
            Index_corr_surf = [Index_corr_surf; Ics+size(V_steelcorr_surf,1)];
            Index_corr_outer = [Index_corr_outer; (1:rebar_point_count{i,1})'+size(V_steelcorr_surf,1)];
            Index_corr_inner = [Index_corr_inner; ((rebar_point_count{i,1}+1):size(Vsc,1))+size(V_steelcorr_surf,1)];
            V_steelcorr_surf = [V_steelcorr_surf; Vsc];
            
            count_steel_points = count_steel_points + steel_point_count{i,1};
            count_steel_F = count_steel_F + steel_F_count{i,1};
            count_corr_points = count_corr_points + rebar_point_count{i,1};
            count_corr_F = count_corr_F + rebar_F_count{i,1};
        end
