function [bar_to_remove, theta, z,ind] = Crack3D(num_longRebar,rebarR,No_ele_rebar_surf,No_ele_rebar_vert,...
    deltax,deltay,deltaz,rebarCenter,exp_z_corr_layer,exp_theta_corr_layer,...
    delta_T,l_A,V_cp_min,V_cp_max,CAR_0,CAR_max,coeff_por,n_non,...
    exp_corr_layer,non_uni_corr_layer,radius,elemlength,conc_no_spacing_bound,...
    X_origin,Y_origin,Z_origin,No_nodes_x,No_nodes_y,No_nodes_z,sub_cell_size,...
    No_nodes_x_rebar,No_nodes_y_rebar,No_nodes_z_rebar,fit_1,E_conc,tensile_strength,...
    E_steel,NU_steel,E_corr,NU_corr,Load_steps,tol_iter,max_iter,flag_solv,DISP,...
    bar_to_remove)

[Points_corr_surf,Points_steel_surf,F_corr_surf,F_steel_surf,ind_rad_steel,...
    rebar_point_count,steel_point_count,rebar_F_count,steel_F_count,...
    top,bottom] = formSteelCorrMesh(num_longRebar,rebarR,No_ele_rebar_surf,...
    No_ele_rebar_vert,deltax,deltay,deltaz,rebarCenter);

[Corr_layer_thick_max,F_diff_non_tot,rebarCenterVec,non_uni_corr_layer] = findCorrDim(num_longRebar,...
    Points_corr_surf,exp_z_corr_layer,exp_theta_corr_layer,rebar_point_count,...
    rebarCenter,delta_T,l_A,rebarR,1,V_cp_min,V_cp_max,CAR_0,CAR_max,...
    coeff_por,n_non,exp_corr_layer,non_uni_corr_layer);

non_uni_corr_layer_comb = [];
for i = 1:num_longRebar
    non_uni_corr_layer_comb = [non_uni_corr_layer_comb; non_uni_corr_layer{i,1}];
end

% MZ move points to match corrosion pattern
[point_steel_surf,corr2steel] = moveCorrPoint(Points_corr_surf,Points_steel_surf,...
    ind_rad_steel,rebarCenterVec,non_uni_corr_layer_comb);

% MZ create surface for concrete meshing
[F_steelcorr_surf,F_corr_surf_final,Index_corr_surf,Index_corr_outer,...
            V_steelcorr_surf] = createRebarSurf(corr2steel,rebarCenterVec,num_longRebar,...
            point_steel_surf,steel_point_count,F_steel_surf,steel_F_count,Points_corr_surf,...
            rebar_point_count,F_corr_surf,rebar_F_count,top,bottom,ind_rad_steel);


% MZ create box for outer concrete boundary
% disp('patching concrete matrix surface')
[F_conc_outer,Points_conc_outer] = patchConcCyl(radius,elemlength,conc_no_spacing_bound);

% MZ generate random points within concrete matrix
% disp('generating random points')
Points_rand = genRandConcCyl(X_origin,Y_origin,Z_origin,radius,Points_conc_outer,...
    elemlength,No_nodes_x,No_nodes_y,No_nodes_z,sub_cell_size,Points_corr_surf,num_longRebar,rebarR);

% MZ alter outer conc surface and combine points and surfaces
[F_conc,Points_conc] = combineConc(V_steelcorr_surf,Points_corr_surf,Points_rand,...
    F_conc_outer,Points_conc_outer,elemlength,radius,rebarR,F_corr_surf,conc_no_spacing_bound);
save('test')
% MZ Delaunay triangulation
[DT_con]=constrainedDelaunayTetGen(Points_conc,F_conc);
TET_con = DT_con.ConnectivityList;
DT_Points_conc = DT_con.Points;
DT_Points_conc = [V_steelcorr_surf; DT_Points_conc((size(Points_corr_surf,1)+1):end,:)];

TET_con(TET_con(:,1) > size(Points_corr_surf,1),1) =...
    TET_con(TET_con(:,1) > size(Points_corr_surf,1),1)...
    + size(V_steelcorr_surf,1) - size(Points_corr_surf,1);
TET_con(TET_con(:,2) > size(Points_corr_surf,1),2) =...
    TET_con(TET_con(:,2) > size(Points_corr_surf,1),2)...
    + size(V_steelcorr_surf,1) - size(Points_corr_surf,1);
TET_con(TET_con(:,3) > size(Points_corr_surf,1),3) =...
    TET_con(TET_con(:,3) > size(Points_corr_surf,1),3)...
    + size(V_steelcorr_surf,1) - size(Points_corr_surf,1);
TET_con(TET_con(:,4) > size(Points_corr_surf,1),4) =...
    TET_con(TET_con(:,4) > size(Points_corr_surf,1),4)...
    + size(V_steelcorr_surf,1) - size(Points_corr_surf,1);

Nodes = DT_Points_conc;

% MZ Construct bar elements
% disp('Constructing bar elements')
Con_bar = constructBar(TET_con);

% disp('Create Corrosion Layer Mesh')
%MZ create mesh of tets for corrosion surf
[DT_corr] = constrainedDelaunayTetGen(V_steelcorr_surf,F_corr_surf_final);

%find new points
new_points = ~ismember(DT_corr.Points,V_steelcorr_surf,'rows');

%determine # new points and add to index
num_new = sum(new_points);
if num_new > 0
    Index_corr_surf_new = zeros(length(Index_corr_surf) + num_new,1);
    Index_corr_surf_new(~new_points,1) = Index_corr_surf;
    Index_corr_surf_new(new_points,1) = (1:num_new)' + length(Nodes);
    Index_corr_surf = Index_corr_surf_new;
    %add in nodes
    Nodes = [Nodes; DT_corr.Points(new_points, :)];
end
Con_quad = Index_corr_surf(DT_corr.ConnectivityList);

% conn_list = DT_corr.ConnectivityList;
% max_c = length(Index_corr_surf);
% conn_list = conn_list(conn_list(:,1)<=max_c & conn_list(:,2)<=max_c & conn_list(:,3)<=max_c & conn_list(:,4)<=max_c,:);
% 
% Con_quad = Index_corr_surf(conn_list);

% disp('Assign corrosion depth to corrosion layer elements')
F_diff_non_elem = modFDiff(Nodes,Index_corr_outer,Con_quad,...
    F_diff_non_tot);

[indexRebar,Nodes] = generateRebarMesh(rebarCenter,Corr_layer_thick_max,...
    rebarR,No_nodes_x_rebar,No_nodes_y_rebar,deltaz,No_nodes_z_rebar,...
    sub_cell_size,point_steel_surf,steel_point_count,F_steel_surf,...
    steel_F_count,Nodes,rebar_point_count,num_longRebar);

% disp("Organizing indices of degrees of freedom")
[~,Ig_bar,Jg_bar] = organizeIndBar(Con_bar);

[~,Ig_tri,Jg_tri] = organizeIndTri(indexRebar);

[~,Ig_quad,Jg_quad] = organizeIndQuad(Con_quad);

% local element stiffness of bar elements
% disp("Construct bar local stiffness matrix")
[K_bar,A_bar,thetaX,thetaY,thetaZ,L_bar] = constructKBar(Nodes,Con_bar,fit_1,E_conc,Ig_bar,Jg_bar,deltaz,No_nodes_z);

tensile_strength_vec = ones(length(A_bar),1) * tensile_strength;
x = [Nodes(Con_bar(:,1),1); Nodes(Con_bar(:,2),1)];
y = [Nodes(Con_bar(:,1),2); Nodes(Con_bar(:,2),2)];

x_cent = mean(x,2) - radius;
y_cent = mean(y,2) - radius;

conc_rad = sqrt(x_cent.^2 + y_cent.^2);

thresh = -1;
s_red = 0;
tensile_strength_vec(conc_rad < (thresh + rebarR),1) = tensile_strength * (1-s_red);

% local stiffness of triangular elements
%MZ define nodes of tet elements
% disp("Construct tet local stiffness matrix")
[K_tri, ~, ~, ~] = constructKTet(Nodes,indexRebar,E_steel,NU_steel,Ig_tri,Jg_tri);

% row 1-48, 48 values for local elements k matrix; column for
% different elements
% local stiffness of quadrilateral elements
% disp("Construct hex local stiffness matrix and force vector")
[K_quad, B, D, V_tet] = constructKTet(Nodes,Con_quad,E_corr,NU_corr,Ig_quad,Jg_quad);
[force_vec, therm_load] = constructFTet(Con_quad,F_diff_non_elem,Nodes,E_corr,NU_corr);

% disp("Constructing displacement vector")
[partion_1,partion_2] = setBound(Nodes);

% partioning of global stiffness matrix
K = K_bar+K_tri+K_quad;

% save('applyLoadStepInput')
[flag_solv,DISP,sigma_bar,bar_to_remove,K,C,res_disp_fin, ind, theta, z] = applyLoadStepMulti(A_bar,...
    E_conc,Load_steps,force_vec,K,partion_1,tol_iter,max_iter,1,...
    Nodes,Con_bar,flag_solv,DISP,thetaX,thetaY,thetaZ,partion_2,...
    tensile_strength,bar_to_remove,Ig_bar,Jg_bar,L_bar,K_tri,K_quad,radius);