function out = MisThresh(grains,axis1,axis2,thresh)

%find the misorientation values from the first 3 areas (must be real data)
delta = GrainsMisori2(grains,axis1);
delta(:,2) = GrainsMisori2(grains,axis2);

[n,~] = size(delta);

for a = 1:n
    if delta(a,2) <= thresh(1,1) || delta(a,2) >= thresh(1,2)
        delta(a,:) = NaN;
    end
end

delta = rmmissing(delta);
out = delta(:,1);