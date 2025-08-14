function Twins = FindTwins(grains,maxT,err)
%identify all grain boundary segments from inputs. "grains" should be a
%grains2d variable, maxT is the maximum twin multiple to be considered, and
%err is the range around the true misorientation to be considered that twin multiple.

%find orientation of the b-axis for each pixel
bori = grains.meanOrientation*Miller(0,0,1,grains.CS);
cori = grains.meanOrientation*Miller(1,0,0,grains.CS);

%count the number of grain boundary pixels between aragonite pixels
gbs = grains.boundary.grainId;
[n_gb,~] = size(gbs);

%check for real boundaries
m_check = isnan(bori.x);

%input cartesian orientation vectors for each grain
b_uv = [bori.x bori.y bori.z];
%calculate the magnitude of the orientation vectors from the first available
%set of numbers (so far they should all be the same magnitude)
for t = 1:n_gb
    if m_check(t,1) == 0
        bmag = sqrt(b_uv(t,1)^2 + b_uv(t,2)^2 + b_uv(t,3)^2);
        break
    end
end
%convert to unit vectors
b_uv = b_uv/bmag;

%input cartesian orientation vectors for each grain
c_uv = [cori.x cori.y cori.z];
%calculate the magnitude of the orientation vectors from the first available
%set of numbers (so far they should all be the same magnitude)
for t = 1:n_gb
    if m_check(t,1) == 0
        cmag = sqrt(c_uv(t,1)^2 + c_uv(t,2)^2 + c_uv(t,3)^2);
        break
    end
end
%convert to unit vectors
c_uv = c_uv/cmag;

%create a matrix to hold the misorientation values, including the IDs of
%the grains with that boundary
deltaB = NaN(n_gb,1);
deltaC = NaN(n_gb,1);

%calculate the misorientation between b-axes of adjacent grains
for q = 1:n_gb
    if gbs(q,1) ~= 0 && gbs(q,2) ~= 0
        deltaB(q,1) = acosd(dot(b_uv(gbs(q,1),:),b_uv(gbs(q,2),:)));
        deltaC(q,1) = acosd(dot(c_uv(gbs(q,1),:),c_uv(gbs(q,2),:)));
    end
end

%ensure all values are =<90° to account for the symmetry of the crystal
for b = 1:n_gb
    if deltaB(b,1) > 90
        deltaB(b,1) = 180 - deltaB(b,1);
    end
    if deltaC(b,1) > 90
        deltaC(b,1) = 180 - deltaB(b,1);
    end
end

Twins = NaN(n_gb,1);

for c = 1:n_gb
    if deltaC(c,1) <= 1
        if deltaB(c,1) > 63.8-err && deltaB(c,1) < 63.8+err
            Twins(c,1) = 1;
        elseif maxT >= 2 && deltaB(c,1) > 52.4-err && deltaB(c,1) < 52.4+err
            Twins(c,1) = 2;
        elseif maxT >= 3 && deltaB(c,1) > 11.4-err && deltaB(c,1) < 11.4+err
            Twins(c,1) = 3;
        elseif maxT >= 4 && deltaB(c,1) > 75.2-err && deltaB(c,1) < 75.2+err
            Twins(c,1) = 4;
        elseif maxT >= 5 && deltaB(c) > 41-err && deltaB(c) < 41+err
            Twins(c) = 5;
        elseif maxT >= 6 && deltaB(c) > 22.8-err && deltaB(c) < 22.8+err
            Twins(c) = 6;
        elseif maxT >= 7 && deltaB(c) > 86.6-err && deltaB(c) < 86.6+err
            Twins(c) = 7;
        elseif maxT >= 8 && deltaB(c) > 29.6-err && deltaB(c) < 29.6+err
            Twins(c) = 8;
        elseif maxT >= 9 && deltaB(c) > 34.2-err && deltaB(c) < 34.2+err
            Twins(c) = 9;
        else
            Twins(c,1) = 0;
        end
    end
end