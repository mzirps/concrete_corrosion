function [me_tri,Ig_tri,Jg_tri] = organizeIndTri(indexRebar)

me_tri = indexRebar';
GetI_tri = @(me_tri) [3*me_tri(1,:)-2; 3*me_tri(1,:)-1; 3*me_tri(1,:);
    3*me_tri(2,:)-2; 3*me_tri(2,:)-1; 3*me_tri(2,:); 3*me_tri(3,:)-2;
    3*me_tri(3,:)-1; 3*me_tri(3,:); 3*me_tri(4,:)-2; 3*me_tri(4,:)-1; 3*me_tri(4,:)];
ii_tri = (1:12)'*ones(1,12);
ii_tri = ii_tri(:);
jj_tri = ones(12,1)*(1:12);
jj_tri = jj_tri(:);

I_tri = GetI_tri(me_tri);

Ig_tri = I_tri(ii_tri,:);
Ig_tri = Ig_tri(:);
Jg_tri = I_tri(jj_tri,:);
Jg_tri = Jg_tri(:);