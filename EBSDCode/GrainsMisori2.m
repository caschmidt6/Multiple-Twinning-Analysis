function out = GrainsMisori2(grains,axis)
%Input "grains" should be a grains2d variable,
%variable "axis" should be a string with a, b, or c, and name should be a
%string with the name of the area. Output the misorientation of the
%selected axis between each pair of adjacent grains.

%convert the axis conventions from the AZtec software to what I'm used to
%working with (the carbonate anions lie in the a-b plane, the b-axis is the
%long axis)
if axis == 'a'
    miller = Miller(0,1,0,grains.CS);
end
if axis == 'b'
    miller = Miller(0,0,1,grains.CS);
end
if axis == 'c'
    miller = Miller(1,0,0,grains.CS);
end

%find orientation of the desired axis for each pixel
ori = grains.meanOrientation*miller;

%identify all grain boundaries between real grains
gId = grains('Aragonite').boundary.grainId;
gId(gId(:,1) == 0,:) = [];
[gId,~,~] = unique(gId,'rows');
ngB = size(gId,1);

%input cartesian orientation vectors for each grain
guv = [ori.x ori.y ori.z];

%calculate the magnitude of the orientation vectors from the first available
%set of numbers (so far they should all be the same magnitude)
for t = 1:ngB
    if isfinite(guv(t,1))
        mag = sqrt(guv(t,1)^2 + guv(t,2)^2 + guv(t,3)^2);
        break
    end
end

%convert to unit vectors
guv = guv/mag;

%create a matrix to hold the misorientation values, including the IDs of
%the grains with that boundary
delta = NaN(ngB,1);

%calculate the misorientation between adjacent grains
for q = 1:ngB
    delta(q,1) = acosd(dot(guv(gId(q,1),:),guv(gId(q,2),:)));
    if delta(q,1) > 90
        delta(q,1) = 180 - delta(q,1);
    end
end

out = delta;