function [me_bar,Ig_bar,Jg_bar] = organizeIndBar(Con_bar)

% indexing of bar elements %MZ concrete lattice element (bar element, 6 dof)
me_bar = Con_bar';
GetI_bar = @(me_bar) [3*me_bar(1,:)-2; 3*me_bar(1,:)-1; 3*me_bar(1,:);
    3*me_bar(2,:)-2; 3*me_bar(2,:)-1; 3*me_bar(2,:)];
ii_bar = (1:6)'*ones(1,6);
ii_bar = ii_bar(:);
jj_bar = ones(6,1)*(1:6);
jj_bar = jj_bar(:);

I_bar = GetI_bar(me_bar);

Ig_bar = I_bar(ii_bar,:);
Ig_bar = Ig_bar(:);
Jg_bar = I_bar(jj_bar,:);
Jg_bar = Jg_bar(:);