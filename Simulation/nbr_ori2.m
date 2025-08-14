function [nbr_cori,nbr_bori] = nbr_ori2(cori,bori,x,y)
%input the orientation matrices and the coordinate you want to check, get
%out a list of neighboring orientations and their weights.

%define size
[L,W] = size(cori);

%define a matrix, nbr_cori, with the orientation of each of the 8
%neighboring pixels as: [+x +x+y +y -x+y -x -x-y -y +x-y]
nbr_cori = nan([1 8]);
nbr_bori = nan([1 8]);
if x ~= L
    nbr_cori(1,1) = cori(x+1,y);
    nbr_bori(1,1) = bori(x+1,y);
    if y ~= W
        nbr_cori(1,2) = cori(x+1,y+1);
        nbr_bori(1,2) = bori(x+1,y+1);
    end
    if y ~= 1
        nbr_cori(1,8) = cori(x+1,y-1);
        nbr_bori(1,8) = bori(x+1,y-1);
    end
end
if x ~= 1
    if y ~= W
        nbr_cori(1,4) = cori(x-1,y+1);
        nbr_bori(1,4) = bori(x-1,y+1);
    end
    nbr_cori(1,5) = cori(x-1,y);
    nbr_bori(1,5) = bori(x-1,y);
    if y ~= 1
        nbr_cori(1,6) = cori(x-1,y-1);
        nbr_bori(1,6) = bori(x-1,y-1);
    end
end
if y ~= W
    nbr_cori(1,3) = cori(x,y+1);
    nbr_bori(1,3) = bori(x,y+1);
end
if y ~= 1
    nbr_cori(1,7) = cori(x,y-1);
    nbr_bori(1,7) = bori(x,y-1);
end