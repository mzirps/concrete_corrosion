function F_diff_non_elem = modFDiff(V_steelcorr_surf,Index_corr_outer,Con_quad,...
            F_diff_non_tot)
        
        %MZ detemine strain based on mesh
        %MZ create indicator vector to determine where 3pt elem vs 1pt elem
        %are located
        indicator = zeros(size(V_steelcorr_surf,1),1);
        indicator(Index_corr_outer,1) = 1;
        indicator = mean(indicator(Con_quad),2);
        
        %MZ extend F_diff_non_tot to be the same length as V_steelcorr_surf
        F_diff_non_tot_ex = zeros(size(V_steelcorr_surf,1),1);
        F_diff_non_tot_ex(Index_corr_outer,1) = F_diff_non_tot;
        
        F_diff_non_elem = zeros(size(Con_quad,1),1);
        
        %MZ create F_diff_non_elem
        F_diff_non_elem(indicator == .75,1) = sum(F_diff_non_tot_ex(Con_quad(indicator == .75,:)),2)/3;
        F_diff_non_elem(indicator == .5,1) = sum(F_diff_non_tot_ex(Con_quad(indicator == .5,:)),2)/2;
        F_diff_non_elem(indicator == .25,1) = sum(F_diff_non_tot_ex(Con_quad(indicator == .25,:)),2);