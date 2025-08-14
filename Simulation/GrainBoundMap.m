function tim = GrainBoundMap(bori,cori,gn,thresh)
%input the simulation outputs and a size threshold to get out a twin
%boundary map made based on the used grain boundaries (excluding those
%below the threshold).

%first remove small grains from consideration
for p = 1:max(gn,[],'all')
    if sum(gn == p,'all') < thresh
        gn(gn == p) = nan;
    end
end

%define size
[L,W] = size(cori);

tw_map = zeros([L W]);

%create the matrix for the b-axis map, then translate the angle into
%grayscale values between 0 and 1 with 0 (black) meaning the b-axis is
%perpendicular to the image plane and 1 (white) the b-axis is parallel.
tim(:,:,1) = sind(bori).^2;
tim(:,:,2) = sind(bori).^2;
tim(:,:,3) = sind(bori).^2;

%loop through every position in the slice and assign the highest adjacent
%twin multiple
for x = 1:L
    for y = 1:W
        if x < L && cori(x,y) == cori(x+1,y) && ~isnan(gn(x,y)) && ~isnan(gn(x+1,y))
            bpt = abs(bori(x,y) - bori(x+1,y));
            if bpt > 90
                bpt = 180 - bpt;
            end
        else
            bpt = 0;
        end
        if y < W && cori(x,y) == cori(x,y+1) && ~isnan(gn(x,y)) && ~isnan(gn(x,y+1))
            cpt = abs(bori(x,y) - bori(x,y+1));
            if cpt > 90
                cpt = 180 - cpt;
            end
        else
            cpt = 0;
        end

        tw_map(x,y) = max([bpt cpt]);

        if tw_map(x,y) ~= 0
            tim(x,y,1:2) = 1;
            tim(x,y,3) = 0;
            if tw_map(x,y) == 64
                tim(x,y,2) = 0;
            elseif tw_map(x,y) == 52
                tim(x,y,2) = 0.5;
            elseif tw_map(x,y) == 76
                tim(x,y,1) = 0;
            end
        end
        if isnan(bori(x,y))
            tim(x,y,:) = 0;
        end
    end
end
