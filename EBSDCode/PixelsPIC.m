function pim = PixelsPIC(ebsd,axis,step)
%Create a single axis orientation map with the PIC mapping color scheme.
%Inputs: "ebsd" must be an EBSD type variable
%"axis" must be a string with the desired axis, 'a', 'b', or 'c'
%"step" should be the step size in µm

%Output: "pim" contains the matrix of RGB values for each pixel, each pixel
%forced to be square.

%convert the axis conventions from the AZtec software to what I'm used to
%working with (the carbonate anions lie in the a-b plane, the b-axis is the
%long axis)
if axis == 'a'
    miller = Miller(0,1,0,ebsd.CS);
end
if axis == 'b'
    miller = Miller(0,0,1,ebsd.CS);
end
if axis == 'c'
    miller = Miller(1,0,0,ebsd.CS);
end

%find orientation of the desired axis for each pixel
ori = ebsd.orientations*miller;

%count the pixels
[n,~] = size(ori);

%convert the in-plane angle (rho) and the out-of-plane angle (theta) into
%hues and values (saturation stays constant at 1), then create matrix to
%house the HSV values of each grain, then fix rounding error
pixel_hue = -1*ori.rho/pi;
pixel_value = ori.theta/pi;
for d = 1:n
    if pixel_value(d,1) > 0.5
        pixel_value(d,1) = 1 - pixel_value(d,1);
    end
end
pixel_value = pixel_value*2;
pixel_hsv = zeros(n,3);
m_check = isnan(pixel_hue);

%fill in hsv values, using 0 instead of NaN, and making the in-plane angle
%symmetrical
for a = 1:n
    if m_check(a,1) == 0
        pixel_hsv(a,2) = 1;
        pixel_hsv(a,3) = pixel_value(a,1);
        if pixel_hue(a,1) < 0
            pixel_hsv(a,1) = -1*pixel_hue(a,1);
        else
            pixel_hsv(a,1) = 1 - pixel_hue(a,1);
        end
    else
        pixel_hsv(a,2) = 0;
        pixel_hsv(a,3) = 1;
    end
end

%convert hsv to rgb colors
pixel_rgb = hsv2rgb(pixel_hsv);

%find the dimensions of the map in pixels
xmax = (max(ebsd.prop.x)/step)+1;
ymax = (max(abs(ebsd.prop.y))/step)+1;

%matrix with the coordinates of each pixel
coords = (ebsd.prop.x/step)+1;
coords(:,2) = abs(ebsd.prop.y/step)+1;

%matrix to hold RGB values per pixel
pim = ones([ymax xmax 3]);

for b = 1:n
    pim(coords(b,2),coords(b,1),1) = pixel_rgb(b,1);
    pim(coords(b,2),coords(b,1),2) = pixel_rgb(b,2);
    pim(coords(b,2),coords(b,1),3) = pixel_rgb(b,3);
end

imshow(pim)