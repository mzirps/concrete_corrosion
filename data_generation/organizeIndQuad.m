function [me_quad,Ig_quad,Jg_quad] = organizeIndQuad(Con_quad)

% indexing of quadrilateral elements %MZ corrosion (24 dof)
me_quad = Con_quad';
GetI_quad = @(me_quad) [3*me_quad(1,:)-2; 3*me_quad(1,:)-1; 3*me_quad(1,:);
    3*me_quad(2,:)-2; 3*me_quad(2,:)-1; 3*me_quad(2,:);
    3*me_quad(3,:)-2; 3*me_quad(3,:)-1; 3*me_quad(3,:);
    3*me_quad(4,:)-2; 3*me_quad(4,:)-1; 3*me_quad(4,:)];
ii_quad = (1:12)'*ones(1,12);
ii_quad = ii_quad(:);
jj_quad = ones(12,1)*(1:12);
jj_quad = jj_quad(:);

I_quad = GetI_quad(me_quad);

Ig_quad = I_quad(ii_quad,:);
Ig_quad = Ig_quad(:);
Jg_quad = I_quad(jj_quad,:);
Jg_quad = Jg_quad(:);