function pic = cori2pic(cori)

%define size
[L,W] = size(cori);

%matrix of hsv values
px_hsv = zeros([L W 3]);

%define saturation and value as 1, as all angles are in-plane
px_hsv(:,:,2) = 1;
px_hsv(:,:,3) = 1;

%convert c-axis orientation values to be between 0 and 1 rather than 0 and
%180, replace nan with 0
px_hsv(:,:,1) = 1-(cori/180);
px_hsv(isnan(cori.*ones(size(px_hsv)))) = 0;


pic = hsv2rgb(px_hsv);
