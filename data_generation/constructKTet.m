function K_tri = constructKTet(Nodes,indexRebar,E_steel,NU_steel,Ig_tri,Jg_tri)

x_1 = Nodes(indexRebar(:,1),1);
x_2 = Nodes(indexRebar(:,2),1);
x_3 = Nodes(indexRebar(:,3),1);
x_4 = Nodes(indexRebar(:,4),1);

y_1 = Nodes(indexRebar(:,1),2);
y_2 = Nodes(indexRebar(:,2),2);
y_3 = Nodes(indexRebar(:,3),2);
y_4 = Nodes(indexRebar(:,4),2);

z_1 = Nodes(indexRebar(:,1),3);
z_2 = Nodes(indexRebar(:,2),3);
z_3 = Nodes(indexRebar(:,3),3);
z_4 = Nodes(indexRebar(:,4),3);

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

D = E_steel/((1+NU_steel)*(1-2*NU_steel))*[1-NU_steel NU_steel NU_steel 0 0 0;
    NU_steel 1-NU_steel NU_steel 0 0 0; NU_steel NU_steel 1-NU_steel 0 0 0;
    0 0 0 (1-2*NU_steel)/2 0 0; 0 0 0 0 (1-2*NU_steel)/2 0; 0 0 0 0 0 (1-2*NU_steel)/2];
D_tri = cell(size(indexRebar,1),1);
D_tri(:) = {D};

% row 1-36, 36 values for local elements k matrix; column for
% different elements
Kg_tri_test = cellfun(@(x,y,z) x'*y*x*z, B, D_tri, num2cell(V_tet), 'UniformOutput', false);
Kg_tri_test = cellfun(@(x) reshape(x, [], 1), Kg_tri_test, 'UniformOutput', false);
Kg_tri_test = cat(1, Kg_tri_test{:});

% global stiffness matrix for tringular elements
nq = size(Nodes,1);
K_tri = sparse(Ig_tri(:),Jg_tri(:),Kg_tri_test(:),3*nq,3*nq);
