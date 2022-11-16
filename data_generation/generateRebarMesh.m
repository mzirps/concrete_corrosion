function [indexRebar,Nodes] = generateRebarMesh(rebarCenter,Corr_layer_thick_max,...
            rebarR,No_nodes_x_rebar,No_nodes_y_rebar,deltaz,No_nodes_z_rebar,...
            sub_cell_size,point_steel_surf,steel_point_count,F_steel_surf,...
            steel_F_count,Nodes,rebar_point_count,num_longRebar)
        
        indexRebarPoint = 1;
        indexRebarIndex = 1;
        countCorrPoint = 1;
        indexRebar  = [];
       
        for i_rebar=1:num_longRebar 
            disp("Generating points and meshing tets for rebar " + string(i_rebar));
            %wu non-corroded steel cell 
            Cell_x_rebar = linspace(rebarCenter(i_rebar,1)-rebarR(i_rebar,1)+Corr_layer_thick_max{i_rebar},...
                rebarCenter(i_rebar,1)+rebarR(i_rebar,1)-Corr_layer_thick_max{i_rebar}, No_nodes_x_rebar+1)';
            Cell_y_rebar = linspace(rebarCenter(i_rebar,2)-rebarR(i_rebar,1)+Corr_layer_thick_max{i_rebar},...
                rebarCenter(i_rebar,2)+rebarR(i_rebar,1)-Corr_layer_thick_max{i_rebar}, No_nodes_y_rebar+1)';
            Cell_z_rebar = linspace(rebarCenter(i_rebar,3), rebarCenter(i_rebar,3) + deltaz(i_rebar,1), No_nodes_z_rebar+1);
            
            % MZ generate random points within steel rebar
            [X_rebar_grid, Y_rebar_grid, Z_rebar_grid] = ndgrid(Cell_x_rebar, Cell_y_rebar, Cell_z_rebar);
            size_x = (2*rebarR(i_rebar,1)- 2*Corr_layer_thick_max{i_rebar})/No_nodes_x_rebar;
            size_y = (2*rebarR(i_rebar,1)- 2*Corr_layer_thick_max{i_rebar})/No_nodes_y_rebar;
            size_z = deltaz(i_rebar,1)/No_nodes_z_rebar;
            X_rand_rebar = (X_rebar_grid + size_x * (1 - sub_cell_size)) + ...
                ((X_rebar_grid + size_x * sub_cell_size) - ...
                (X_rebar_grid + size_x * (1 - sub_cell_size))) .* rand(size(X_rebar_grid,1),size(X_rebar_grid,2),size(X_rebar_grid,3));
            Y_rand_rebar = (Y_rebar_grid + size_y * (1 - sub_cell_size)) + ...
                ((Y_rebar_grid + size_y * sub_cell_size) - ...
                (Y_rebar_grid + size_y * (1 - sub_cell_size))) .* rand(size(Y_rebar_grid,1),size(Y_rebar_grid,2),size(Y_rebar_grid,3));
            Z_rand_rebar = (Z_rebar_grid + size_z * (1 - sub_cell_size)) + ...
                ((Z_rebar_grid + size_z * sub_cell_size) - ...
                (Z_rebar_grid + size_z * (1 - sub_cell_size))) .* rand(size(Z_rebar_grid,1),size(Z_rebar_grid,2),size(Z_rebar_grid,3));
            X_rand_rebar = X_rand_rebar(1:end-1,1:end-1,1:end-1);
            Y_rand_rebar = Y_rand_rebar(1:end-1,1:end-1,1:end-1);
            Z_rand_rebar = Z_rand_rebar(1:end-1,1:end-1,1:end-1);
            Points_rand_rebar = [X_rand_rebar(:) Y_rand_rebar(:) Z_rand_rebar(:)];
        
            inCircleRebar = logical(1 - checkInCircle3D(Points_rand_rebar, point_steel_surf(indexRebarPoint:(steel_point_count{i_rebar,1}+(indexRebarPoint-1)),:), num_longRebar));
            Points_rand_rebar = Points_rand_rebar(inCircleRebar,:);
            
            Points = [point_steel_surf(indexRebarPoint:(steel_point_count{i_rebar,1}+(indexRebarPoint-1)),:); Points_rand_rebar];
            Index = F_steel_surf(indexRebarIndex:(steel_F_count{i_rebar,1}+(indexRebarIndex-1)),:)-(indexRebarPoint-1);
            [DT_rebar]=constrainedDelaunayTetGen(Points,Index);
            
            indexRow1 = zeros(size(DT_rebar.ConnectivityList,1), 1);
            indexRow2 = zeros(size(DT_rebar.ConnectivityList,1), 1);
            indexRow3 = zeros(size(DT_rebar.ConnectivityList,1), 1);
            indexRow4 = zeros(size(DT_rebar.ConnectivityList,1), 1);
            
            connectRow1 = DT_rebar.ConnectivityList(:,1);
            connectRow2 = DT_rebar.ConnectivityList(:,2);
            connectRow3 = DT_rebar.ConnectivityList(:,3);
            connectRow4 = DT_rebar.ConnectivityList(:,4);
            
            
            %MZ adjust indexing for inner  points
            indexRow1(connectRow1 > steel_point_count{i_rebar,1}) = connectRow1(connectRow1 > steel_point_count{i_rebar,1}, 1) + size(Nodes,1) - steel_point_count{i_rebar,1};
            indexRow2(connectRow2 > steel_point_count{i_rebar,1}) = connectRow2(connectRow2 > steel_point_count{i_rebar,1}, 1) + size(Nodes,1) - steel_point_count{i_rebar,1};
            indexRow3(connectRow3 > steel_point_count{i_rebar,1}) = connectRow3(connectRow3 > steel_point_count{i_rebar,1}, 1) + size(Nodes,1) - steel_point_count{i_rebar,1};
            indexRow4(connectRow4 > steel_point_count{i_rebar,1}) = connectRow4(connectRow4 > steel_point_count{i_rebar,1}, 1) + size(Nodes,1) - steel_point_count{i_rebar,1};
            %MZ adjust indexing for surface points
            indexRow1(connectRow1 <= steel_point_count{i_rebar,1}) = connectRow1(connectRow1 <= steel_point_count{i_rebar,1},1) + rebar_point_count{i_rebar,1} + (countCorrPoint-1) + (indexRebarPoint-1);
            indexRow2(connectRow2 <= steel_point_count{i_rebar,1}) = connectRow2(connectRow2 <= steel_point_count{i_rebar,1},1) + rebar_point_count{i_rebar,1} + (countCorrPoint-1) + (indexRebarPoint-1);
            indexRow3(connectRow3 <= steel_point_count{i_rebar,1}) = connectRow3(connectRow3 <= steel_point_count{i_rebar,1},1) + rebar_point_count{i_rebar,1} + (countCorrPoint-1) + (indexRebarPoint-1);
            indexRow4(connectRow4 <= steel_point_count{i_rebar,1}) = connectRow4(connectRow4 <= steel_point_count{i_rebar,1},1) + rebar_point_count{i_rebar,1} + (countCorrPoint-1) + (indexRebarPoint-1);
            
            indexRebar = [indexRebar; indexRow1 indexRow2 indexRow3 indexRow4];
            Nodes = [Nodes; DT_rebar.Points(steel_point_count{i_rebar,1}+1:end,:)];
            
            indexRebarPoint = indexRebarPoint + steel_point_count{i_rebar,1};
            indexRebarIndex = indexRebarIndex + steel_F_count{i_rebar,1};
            countCorrPoint = countCorrPoint +  rebar_point_count{i_rebar,1};
            
        end