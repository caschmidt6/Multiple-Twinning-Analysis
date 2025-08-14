function out = MisThresh(grains,axis1,axis2,thresh)
%"grains" should be a grains2d variable, "axis1" and "axis2" should both be
%strings with "a", "b", or "c", and thresh should be a bracketed range (ex.
%[0 1].

%Output the misorientation of axis1 between each pair of adjacent grains
%that falls within the given "thresh" range of misorientation of axis2.

%find the misorientation values of each pair of adjacent grains for both
%axes
delta = GrainsMisori2(grains,axis1);
delta(:,2) = GrainsMisori2(grains,axis2);

[n,~] = size(delta);

%remove the data for adjacent grains that have an axis2 misorientation
%outside the given range
for a = 1:n
    if delta(a,2) <= thresh(1,1) || delta(a,2) >= thresh(1,2)
        delta(a,:) = NaN;
    end
end

delta = rmmissing(delta);
out = delta(:,1);