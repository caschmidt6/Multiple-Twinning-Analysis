function [hist,hc] = TwinBHist(cori,bori,gn,thresh)
%input the b- and c- orientation maps and the grain ID matrix to output
%misorientation between adjacent grains co-oriented in c-axis orientation.
%cori and bori are axial orientation matrices, gn is the grain ID matrix,
%and thresh is the minimum area of pixels a grain must have to be counted.

%first remove small grains from consideration
for p = 1:max(gn,[],'all')
    if sum(gn == p,'all') < thresh
        gn(gn == p) = nan;
    end
end

%define size of matrix and which pixels are filled
[L,W] = size(cori);
%use gn for the index so that small grains don't get counted
ind = isfinite(gn);

%"dif" will contain all "misorientation" values between adjacent pixels not
%in the same grain but co-oriented in c, the first two columns will
%containt the grains being compared and the third will contain the
%misorientation value. Grain IDs always sorted smallest to largest to avoid
%comparing the same grain twice.
dif = NaN(L*W*2,3);
%counter for loop
m = 1;

%loop through each pixel, and determine if it and its positive direction neighbors
%1. exist 2. are in different grains and 3. are co-oriented in c-axis orientation
for b = 1:L
    for c = 1:W
        if ind(b,c) == 1
            %first check neighbor in the +x direction
            if b < L && ind(b+1,c) == 1 && gn(b,c) ~= gn(b+1,c) && cori(b,c) == cori(b+1,c)
                dif(m,1) = min([gn(b,c) gn(b+1,c)]);
                dif(m,2) = max([gn(b,c) gn(b+1,c)]);
                dif(m,3) = abs(bori(b,c) - bori(b+1,c));
                if dif(m,3) > 90
                    dif(m,3) = 180 - dif(m,3);
                end
                m = m+1;
            end
            %next check neighbor in +y direction
            if c < W && ind(b,c+1) == 1 && gn(b,c) ~= gn(b,c+1) && cori(b,c) == cori(b,c+1)
                dif(m,1) = min([gn(b,c) gn(b,c+1)]);
                dif(m,2) = max([gn(b,c) gn(b,c+1)]);
                dif(m,3) = abs(bori(b,c) - bori(b,c+1));
                if dif(m,3) > 90
                    dif(m,3) = 180 - dif(m,3);
                end
                m = m+1;
            end
        end
    end
end

%remove missing rows
dif = rmmissing(dif);

%find unique rows
dif = unique(dif,'rows');

%output just the misorientations and the twinning frequencies
hist = dif(:,3);
hc = TwinFrequencySim(hist);
