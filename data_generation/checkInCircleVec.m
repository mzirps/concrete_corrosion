function inCircle=checkInCircleVec(x,y,rebarCenterX,rebarCenterY,rebarR)
    inCircle=zeros(size(x,1),1);
    inCircle((sqrt((x-rebarCenterX).^2 + (y-rebarCenterY).^2)) <= rebarR,1) = 1;