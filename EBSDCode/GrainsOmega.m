function out = GrainsOmega(grains)
%input a grains2d variable, get out the average misorientation (omega)
%between each pair of adjacent grains.

%identify all grain boundaries between real grains
gId = grains('Aragonite').boundary.grainId;
gId(gId(:,1) == 0,:) = [];
[gId,~,~] = unique(gId,'rows');
ngB = size(gId,1);

%find orientations of each grain
ori = grains.meanOrientation;

omega = nan([ngB 1]);

for a = 1:ngB
    omega(a) = angle(ori(gId(a,1)),ori(gId(a,2)))./degree;
    if omega(a) > 90
        omega(a) = 180 - omega(a);
    end
end

out = rmmissing(omega);