function force_vec = constructFTet(Con_quad,F_diff_non_elem,Nodes,E_corr,NU_corr)

x_1 = Nodes(Con_quad(:,1),1);
y_1 = Nodes(Con_quad(:,1),2);
z_1 = Nodes(Con_quad(:,1),3);
x_2 = Nodes(Con_quad(:,2),1);
y_2 = Nodes(Con_quad(:,2),2);
z_2 = Nodes(Con_quad(:,2),3);
x_3 = Nodes(Con_quad(:,3),1);
y_3 = Nodes(Con_quad(:,3),2);
z_3 = Nodes(Con_quad(:,3),3);
x_4 = Nodes(Con_quad(:,4),1);
y_4 = Nodes(Con_quad(:,4),2);
z_4 = Nodes(Con_quad(:,4),3);


%MZ determine all alpha, beta, gamma, and delta
beta_1 = -(y_3.*z_4 + y_2.*z_3 + z_2.*y_4 - y_2.*z_4 - z_3.*y_4 -  z_2.*y_3);

gamma_1 = x_3.*z_4 + x_2.*z_3 + z_2.*x_4 - x_2.*z_4 - z_3.*x_4 - z_2.*x_3;

delta_1 =-(x_3.*y_4 + x_2.*y_3 + y_2.*x_4 - x_2.*y_4 - y_3.*x_4 - y_2.*x_3);

beta_2 = y_3.*z_4 + y_1.*z_3 + z_1.*y_4 - y_1.*z_4 - z_3.*y_4 -  z_1.*y_3;

gamma_2 = -(x_3.*z_4 + x_1.*z_3 + z_1.*x_4 - x_1.*z_4 - z_3.*x_4 - z_1.*x_3);

delta_2 = x_3.*y_4 + x_1.*y_3 + y_1.*x_4 - x_1.*y_4 - y_3.*x_4 - y_1.*x_3;

beta_3 = -(y_2.*z_4 + y_1.*z_2 + z_1.*y_4 - y_1.*z_4 - z_2.*y_4 -  z_1.*y_2);

gamma_3 = x_2.*z_4 + x_1.*z_2 + z_1.*x_4 - x_1.*z_4 - z_2.*x_4 - z_1.*x_2;

delta_3 = -(x_2.*y_4 + x_1.*y_2 + y_1.*x_4 - x_1.*y_4 - y_2.*x_4 - y_1.*x_2);

beta_4 = y_2.*z_3 + y_1.*z_2 + z_1.*y_3 - y_1.*z_3 - z_2.*y_3 -  z_1.*y_2;

gamma_4 =-(x_2.*z_3 + x_1.*z_2 + z_1.*x_3 - x_1.*z_3 - z_2.*x_3 - z_1.*x_2);

delta_4 = x_2.*y_3 + x_1.*y_2 + y_1.*x_3 - x_1.*y_3 - y_2.*x_3 - y_1.*x_2;

V_tet = (1/6)*(x_1.*y_3.*z_2 - x_1.*y_2.*z_3 + x_2.*y_1.*z_3...
    - x_2.*y_3.*z_1 - x_3.*y_1.*z_2 + x_3.*y_2.*z_1 + x_1.*y_2.*z_4...
    - x_1.*y_4.*z_2 - x_2.*y_1.*z_4 + x_2.*y_4.*z_1 + x_4.*y_1.*z_2...
    - x_4.*y_2.*z_1 - x_1.*y_3.*z_4 + x_1.*y_4.*z_3 + x_3.*y_1.*z_4...
    - x_3.*y_4.*z_1 - x_4.*y_1.*z_3 + x_4.*y_3.*z_1 + x_2.*y_3.*z_4...
    - x_2.*y_4.*z_3 - x_3.*y_2.*z_4 + x_3.*y_4.*z_2 + x_4.*y_2.*z_3 - x_4.*y_3.*z_2);

B = cellfun(@(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,v)...
    (1/(6*v))*[x1 0 0 x2 0 0 x3 0 0 x4 0 0;
    0 y1 0 0 y2 0 0 y3 0 0 y4 0; 0 0 z1 0 0 z2 0 0 z3 0 0 z4;
    y1 x1 0 y2 x2 0 y3 x3 0 y4 x4 0; 0 z1 y1 0 z2 y2 0 z3 y3 0 z4 y4;
    z1 0 x1 z2 0 x2 z3 0 x3 z4 0 x4], num2cell(beta_1)...
    , num2cell(beta_2), num2cell(beta_3), num2cell(beta_4),...
    num2cell(gamma_1), num2cell(gamma_2), num2cell(gamma_3)...
    , num2cell(gamma_4), num2cell(delta_1), num2cell(delta_2)...
    , num2cell(delta_3), num2cell(delta_4), num2cell(V_tet)...
    , 'UniformOutput', false);

D = E_corr/((1+NU_corr)*(1-2*NU_corr))*[1-NU_corr NU_corr NU_corr 0 0 0;
    NU_corr 1-NU_corr NU_corr 0 0 0; NU_corr NU_corr 1-NU_corr 0 0 0;
    0 0 0 (1-2*NU_corr)/2 0 0; 0 0 0 0 (1-2*NU_corr)/2 0; 0 0 0 0 0 (1-2*NU_corr)/2];
D_hex = cell(size(Con_quad,1),1);
D_hex(:) = {D};

% define force vector
force_vec = zeros(3*size(Nodes,1),1);
eps_therm = cellfun(@(x) [x; x; x; 0; 0; 0], num2cell(F_diff_non_elem),'UniformOutput', false);

therm_load = cellfun(@(a,b,c,d) a'*b*c*d, B, D_hex, eps_therm, num2cell(V_tet), 'UniformOutput',false);

%reorder all force vectors and add them onto appropriate global
%nodes for following two paragraph
force_vec_therm_ind(1,:) = Con_quad(:,1)'*3-2;
force_vec_therm_ind(2,:) = Con_quad(:,1)'*3-1;
force_vec_therm_ind(3,:) = Con_quad(:,1)'*3;
force_vec_therm_ind(4,:) = Con_quad(:,2)'*3-2;
force_vec_therm_ind(5,:) = Con_quad(:,2)'*3-1;
force_vec_therm_ind(6,:) = Con_quad(:,2)'*3;
force_vec_therm_ind(7,:) = Con_quad(:,3)'*3-2;
force_vec_therm_ind(8,:) = Con_quad(:,3)'*3-1;
force_vec_therm_ind(9,:) = Con_quad(:,3)'*3;
force_vec_therm_ind(10,:) = Con_quad(:,4)'*3-2;
force_vec_therm_ind(11,:) = Con_quad(:,4)'*3-1;
force_vec_therm_ind(12,:) = Con_quad(:,4)'*3;

therm_load_vec = cat(1,therm_load{:});
force_vec_therm = [therm_load_vec reshape(force_vec_therm_ind,[],1)];
force_vec_therm_fin = accumarray(force_vec_therm(:,2), force_vec_therm(:,1));
force_vec_ind = (1:size(force_vec_therm_fin,1))';
force_vec(force_vec_ind,1) = force_vec(force_vec_ind,1)+force_vec_therm_fin;