function [nbr_gn] = nbr_gn(gn,x,y)
%input the orientation matrices and the coordinate you want to check, get
%out a list of neighboring orientations and their weights.

%define size
[L,W] = size(gn);

%define a matrix, nbr_cori, with the orientation of each of the 8
%neighboring pixels as: [+x +x+y +y -x+y -x -x-y -y +x-y]
nbr_gn = nan([1 8]);
if x ~= L
    nbr_gn(1,1) = gn(x+1,y);
    if y ~= W
        nbr_gn(1,2) = gn(x+1,y+1);
    end
    if y ~= 1
        nbr_gn(1,8) = gn(x+1,y-1);
    end
end
if x ~= 1
    if y ~= W
        nbr_gn(1,4) = gn(x-1,y+1);
    end
    nbr_gn(1,5) = gn(x-1,y);
    if y ~= 1
        nbr_gn(1,6) = gn(x-1,y-1);
    end
end
if y ~= W
    nbr_gn(1,3) = gn(x,y+1);
end
if y ~= 1
    nbr_gn(1,7) = gn(x,y-1);
end