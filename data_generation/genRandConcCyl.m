function Points_rand = genRandConcCyl(X_origin,Y_origin,Z_origin,radius,Points_conc_outer,...
            elemlength,No_nodes_x,No_nodes_y,No_nodes_z,sub_cell_size,Points_corr_surf,num_longRebar,rebarR)
        
        Cell_x_conc = linspace(X_origin, radius*2, No_nodes_x+1)';
        Cell_y_conc = linspace(Y_origin, radius*2, No_nodes_y+1)';
        Cell_z_conc = linspace(Z_origin, elemlength, No_nodes_z+1)';
        [X_conc_grid, Y_conc_grid, Z_conc_grid] = ndgrid(Cell_x_conc, Cell_y_conc, Cell_z_conc);
        size_x = (radius*2)/No_nodes_x;
        size_y = (radius*2)/No_nodes_y;
        size_z = elemlength/No_nodes_z;
        X_rand_conc = (X_conc_grid + size_x * (1 - sub_cell_size)) + ...
            ((X_conc_grid + size_x * sub_cell_size) - ...
            (X_conc_grid + size_x * (1 - sub_cell_size))) .* rand(size(X_conc_grid,1),size(X_conc_grid,2),size(X_conc_grid,3));
        Y_rand_conc = (Y_conc_grid + size_y * (1 - sub_cell_size)) + ...
            ((Y_conc_grid + size_y * sub_cell_size) - ...
            (Y_conc_grid + size_y * (1 - sub_cell_size))) .* rand(size(Y_conc_grid,1),size(Y_conc_grid,2),size(Y_conc_grid,3));
        Z_rand_conc = (Z_conc_grid + size_z * (1 - sub_cell_size)) + ...
            ((Z_conc_grid + size_z * sub_cell_size) - ...
            (Z_conc_grid + size_z * (1 - sub_cell_size))) .* rand(size(Z_conc_grid,1),size(Z_conc_grid,2),size(Z_conc_grid,3));
        X_rand_conc = X_rand_conc(1:end-1,1:end-1,1:end-1);
        Y_rand_conc = Y_rand_conc(1:end-1,1:end-1,1:end-1);
        Z_rand_conc = Z_rand_conc(1:end-1,1:end-1,1:end-1);
        Points_rand = [X_rand_conc(:) Y_rand_conc(:) Z_rand_conc(:)];
        
        % create larger cylinder to remove points within corrosion surface
        input.cylRadius = rebarR + 10^-4;
        input.pointSpacing = .001;
        input.cylHeight = elemlength+.2;
        [~,closedCorrSurf,~] = patchClosedCylinder(input);
        closedCorrSurf = [closedCorrSurf(:,1:2)+radius closedCorrSurf(:,3) + elemlength/2];
        inCircle = checkInCircle3D(Points_rand, closedCorrSurf, num_longRebar);
        Points_rand = Points_rand(inCircle,:);
        
        % create larger cylinder to remove points within conc surface
        input.cylRadius = radius - 10^-4;
        input.pointSpacing = .001;
        input.cylHeight = elemlength+.2;
        [~,closedConcSurf,~] = patchClosedCylinder(input);
        closedConcSurf = [closedConcSurf(:,1:2)+radius closedConcSurf(:,3) + elemlength/2];
        inConc = logical(1 - checkInCircle3D(Points_rand, closedConcSurf, num_longRebar));
        Points_rand = Points_rand(inConc,:);