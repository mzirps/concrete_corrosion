function Con_bar = constructBar(TET_con)

        Con_bar_1 = [TET_con(:,1) TET_con(:,2); TET_con(:,2) TET_con(:,3);...
            TET_con(:,3) TET_con(:,4); TET_con(:,4) TET_con(:,1)];
        Con_bar_2 = sort(Con_bar_1,2);
        Con_bar_2(:,3) = (1:1:length(Con_bar_2))';
        
        Con_bar_3 = accumarray(Con_bar_2(:,1:2), Con_bar_2(:,3)', [], @prod, 0, true);
        [i, j, ~] = find(Con_bar_3);
        Con_bar = [i, j];