function [F_conc_outer,Points_conc_outer] = patchConcCyl(radius,elemlength,conc_no_spacing_bound)

input.cylRadius = radius;
input.numRadial = round((2*radius*pi())/conc_no_spacing_bound);
input.cylHeight = elemlength;
input.numHeight = round(elemlength/conc_no_spacing_bound);
input.meshType = 'tri';
input.closeOpt = 1;

[F_conc_outer, Points_conc_outer, ~] = patchcylinder(input);
Points_conc_outer = [Points_conc_outer(:,1) + radius Points_conc_outer(:,2) + radius Points_conc_outer(:,3) + elemlength/2];
