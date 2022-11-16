function [inCircle]=checkInCircle3D(Points_rand, Points_rebar_surf, num_longRebar)
    step = size(Points_rebar_surf,1)/num_longRebar;
    startstep = 1;
    totalstep = step;
    inCircle = zeros(size(Points_rand,1),num_longRebar);
    for i = 1:num_longRebar
        rebar_points = Points_rebar_surf(startstep:totalstep,:);
        inCircle(:,i) = inhull(Points_rand, rebar_points);
        startstep = startstep + step;
        totalstep = totalstep + step;
    end
   inCircle = -(inCircle - 1);
   inCircle = all(inCircle,2);
   
   
    