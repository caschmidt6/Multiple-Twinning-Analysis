function ipf = PixelsIPF(ebsd,step)
%input an ebsd variable, output the matrix of RGB values associated with
%the IPF-z map. "step" is the step size, tells the code how big a pixel is.

%define color key for aragonite
key = ipfColorKey(ebsd('aragonite'));

%convert orientations to rgb colors via ipf color key
rgb = key.orientation2color(ebsd.orientations);

%find the dimensions of the map in pixels
xmax = (max(ebsd.prop.x)/step)+1;
ymax = (max(abs(ebsd.prop.y))/step)+1;

%matrix with the coordinates of each pixel
coords = (ebsd.prop.x/step)+1;
coords(:,2) = abs(ebsd.prop.y/step)+1;

%matrix to hold RGB values per pixel
ipf = ones([ymax xmax 3]);

%define number of pixels
n = size(rgb,1);

%loop through rgb matrix and input values into ipf map
for a = 1:n
    ipf(coords(a,2),coords(a,1),:) = rgb(a,:);
end