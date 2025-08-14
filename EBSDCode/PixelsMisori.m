function out = PixelsMisori(ebsd,axis)
%must run EBSDimportscript.m first. Variable "ebsd" remains "ebsd",
%variable "axis" should be a string with a, b, or c, and name should be a
%string with the name of the area.


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

%create a matrix with the actual dimensions of the image taken, which will
%house the x, y, and z parts of each unit vector (uv). Must first find the
%step size, as x and y coordinates are in µm, not pixels. 
stepsize = min(ebsd.prop.x(ebsd.prop.x>0));
xcoord = abs(ebsd.prop.x)/stepsize;
ycoord = abs(ebsd.prop.y)/stepsize;
xmax = max(xcoord) + 1;
ymax = max(ycoord) + 1;
uv_map = zeros(xmax,ymax,3);

%find number of total pixels
[n,~] = size(ori.x);

%check if pixel contains data
m_check = isnan(ori.x);

%populate the unit vector map, replacing NaN with 0
for j = 1:n
    if m_check(j,1) == 0
        uv_map(xcoord(j,1)+1,ycoord(j,1)+1,1) = ori.x(j,1);
        uv_map(xcoord(j,1)+1,ycoord(j,1)+1,2) = ori.y(j,1);
        uv_map(xcoord(j,1)+1,ycoord(j,1)+1,3) = ori.z(j,1);
    end
end

%the oreintation vectors are not unit vectors yet, so find the magnitude
%and divide by it
for t = 1:n
    if m_check(t,1) == 0
        mag = sqrt(ori.x(t,1)^2 + ori.y(t,1)^2 + ori.z(t,1)^2);
        break
    end
end

uv_map = uv_map/mag;

%create matrix to hold misorientation values
deltaC = NaN(n,1);
%counter for upcoming loop
a = 1;

%loop through the map, noting the angle between each pair of pixels
for x = 1:xmax
    for y = 1:ymax
        %ignore un-indexed pixels
        if uv_map(x,y,1) ~= 0 && uv_map(x,y,2) ~= 0 && uv_map(x,y,3) ~= 0
            %measure misorientation with the pixel below
            if x < xmax
                if uv_map(x+1,y,1) ~= 0 && uv_map(x+1,y,2) ~= 0 && uv_map(x+1,y,3) ~= 0
                    %measure the misorientation with the dot product,
                    %(arccos(a•b))
                    deltaC(a,1) = acos(dot([uv_map(x,y,1),uv_map(x,y,2),uv_map(x,y,3)],[uv_map(x+1,y,1),uv_map(x+1,y,2),uv_map(x+1,y,3)]));
                    a = a + 1;
                end
            end
            %measure misorientation with the pixel to the right
            if y < ymax
                if uv_map(x,y+1,1) ~= 0 && uv_map(x,y+1,2) ~= 0 && uv_map(x,y+1,3) ~= 0
                    deltaC(a,1) = acos(dot([uv_map(x,y,1),uv_map(x,y,2),uv_map(x,y,3)],[uv_map(x,y+1,1),uv_map(x,y+1,2),uv_map(x,y+1,3)]));
                    a = a + 1;
                end
            end
        end
    end
end

%remove excess NaN from deltaC
deltaC = rmmissing(deltaC);

%convert deltaC into degrees
deltaC = deltaC*(180/pi);

%ensure all values are =<90° to account for the symmetry of the crystal
a = a-1;
for b = 1:a
    if deltaC(b,1) > 90
        deltaC(b,1) = 180 - deltaC(b,1);
    end
end

out = deltaC;