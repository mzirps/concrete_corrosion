function [K_bar,A_bar,thetaX,thetaY,thetaZ,L_bar] = constructKBar(Nodes,Con_bar,fit_1,E_conc,Ig_bar,Jg_bar)

Delta_x_bar = Nodes(Con_bar(:,2),1)-Nodes(Con_bar(:,1),1);
Delta_y_bar = Nodes(Con_bar(:,2),2)-Nodes(Con_bar(:,1),2);
Delta_z_bar = Nodes(Con_bar(:,2),3)-Nodes(Con_bar(:,1),3);
L_bar = sqrt(Delta_x_bar.^2+Delta_y_bar.^2+Delta_z_bar.^2);

A_bar = (delta_z./No_nodes_z).*fit_1.*L_bar;
thetaX = Delta_x_bar./L_bar;
thetaY = Delta_y_bar./L_bar;
thetaZ = Delta_z_bar./L_bar;

const = E_conc*A_bar./L_bar;
kc_prime = cell(size(Con_bar,1),1);

Tstar = cellfun(@(x,y,z) [x 0; y 0; z 0; 0 x; 0 y; 0 z], num2cell(thetaX), num2cell(thetaY), num2cell(thetaZ), 'UniformOutput', false);
kc_prime(:) = {[1 -1; -1 1]};
Kg_bar_m = cellfun(@(x,y,z) z*x*y*x', Tstar, kc_prime, num2cell(const), 'UniformOutput',false);

Kg_bar = cellfun(@(x) reshape(x, [], 1), Kg_bar_m, 'UniformOutput', false);
Kg_bar = cat(1, Kg_bar{:});

% global stiffness matrix for bar elements
nq = size(Nodes,1);
K_bar = sparse(Ig_bar(:),Jg_bar(:),Kg_bar(:),3*nq,3*nq); %wu construct global matrix
