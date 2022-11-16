function [deltax,deltay,deltaz,rebarR,rebarCenter] = getRebarDim(longRebar,num_longRebar)

totalrebarSet=size(longRebar{2},1);
totallongRebar=1;
deltax = zeros(num_longRebar,1);
deltay = zeros(num_longRebar,1);
deltaz = zeros(num_longRebar,1);
rebarR=zeros(num_longRebar,1);
rebarCenter=zeros(num_longRebar,3);

for k=1:totalrebarSet
    directionx=(longRebar{2}(k,7)-longRebar{2}(k,5));
    directiony=(longRebar{2}(k,8)-longRebar{2}(k,6));
    for j=1:longRebar{2}(k,1)
        deltax(totallongRebar,1) = longRebar{1}(k,4)-longRebar{1}(k,1);
        deltay(totallongRebar,1) = longRebar{1}(k,5)-longRebar{1}(k,2);
        deltaz(totallongRebar,1) = (longRebar{1}(k,6)-longRebar{1}(k,3));
        rebarR(totallongRebar,1)=longRebar{2}(k,2)/2;
        if longRebar{2}(k,1) > 1
            rebarCenter(totallongRebar,:) = [longRebar{1}(k,1)+directionx*(j-1)/(longRebar{2}(k,1)-1)...
                longRebar{1}(k,2)+directiony*(j-1)/(longRebar{2}(k,1)-1) longRebar{1}(totalrebarSet,3)];
        else
            rebarCenter(totallongRebar,:) = [longRebar{1}(k,1) longRebar{1}(k,2) longRebar{1}(totalrebarSet,3)];
        end
        totallongRebar=totallongRebar+1;
    end
end