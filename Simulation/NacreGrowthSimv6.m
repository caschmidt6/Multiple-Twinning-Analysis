function [pim,tim,bori,cori,gn] = NacreGrowthSimv6(NP,GP,TP,L,W,linspc,osz,ospc,ostg,t,nameC,nameA)
%% Inputs
%"NP" is the probability of nucleating a completely new orientation

%"GP" is the probability of growth along the a- or b- axes

%"TP" is the probability of twinning when growth occurs in the ab-plane

%"L" and "W" are the dimensions of the simulation

%"linspc" is the spacing of the simulated organic sheets

%"osz" is the size of the openings in the simulated organic sheets, "ospc"
%is the spacing between openings, and "ostg" is distance each opening is
%staggered from the next sheet (for simulating sheet nacre)

%"t" is the number of time intervals to simulate

%"nameC" and "nameA" are the names of the .mp4 files created with the
%movies of the simulated growth. Names must be strings that end with .mp4.
%Recording the movies increases runtime significantly, so recording can be
%suppressed by inputting 0 (as a number, not a string)
%% Set-up
%matrices to hold orientations.
cori = nan([L W]);
bori = nan([L W]);

%matrix to hold grain numbers (gn) and a counter to keep track of how many
%grains there are
gn = nan([L W]);

%fill the bottom line with random orientations
cori(L,:) = randi(180,[1 W]);
bori(L,:) = randi(180,[1 W]);

%assign each nucleation point a grain ID
gn(L,:) = 1:W;
gc = W;

%matrix to keep track of which positions are filled, 1 is filled, 0 is
%empty
ind = isfinite(cori);

%list filled coordinates
[xfil,yfil] = find(ind);
nf = size(xfil,1);

%matrix to hold probability of growth "pg" for each pixel in matrix,
%which also keeps track of the possible locations for
%growth/nucleation.
pg = zeros([L W]);

%determine positions where growth/nucleation are not allowed, first
%calculate the number of lines between the center line and the edge of the
%simulation
nlin = floor((L-100)/(linspc));

%matrix to store the line placements
clin = zeros([L W]);

%find x values for each line
xlin = linspc*[1:nlin];

%place the lines, spaced from the center
clin(xlin',:) = 1;

%find number of openings per line
no = floor(W/(osz+ospc));
if (W - no*(osz+ospc)) < ospc
    no = no+1;
end

%find the y positions of the openings if you start at 1
ylin = zeros([nlin no*osz]);
spc = 1:osz;
for j = 1:no
    ylin(1,(j-1)*osz+1:j*osz) = spc + (j-1)*(ospc + osz);
end

%find the y positions of the openings for each line moving upward, then
%correct the values to wrap around to the beginning. 
for u = 1:nlin-1
    ylin(u+1,:) = ylin(u,:) + ostg;
    for v = 1:no*osz
        if ylin(u+1,v) > W
            ylin(u+1,v) = ylin(u+1,v) - W;
        end
    end
end

%for each line, set the openings in the matrix
for r = 1:nlin
    for s = 1:no*osz
        clin(xlin(1,r),ylin(r,s)) = 0;
    end
end


%matrix to hold list of neighbors, formatted as an LxLx8 matrix where x
%and y are the coordinates of the pixel, and z are the 8 orientations
%of the neighbors in this order from 1 to 8: 1=+x 2=+x+y 3=+y 4=-x+y 
%5=-x 6=-x-y 7=-y 8=+x-y
nbr_b = nan([L W 8]);
nbr_c = nan([L W 8]);
nbr_g = nan([L W 8]);

%find initial coordinates for addition
for a = 1:nf
    if xfil(a) ~= L
        %check +x
        if ind(xfil(a)+1,yfil(a)) == 0 && clin(xfil(a)+1,yfil(a)) == 0
            pg(xfil(a)+1,yfil(a)) = pg(xfil(a)+1,yfil(a)) + GP*(10*sind(cori(xfil(a,1),yfil(a,1)))^2+cosd(cori(xfil(a,1),yfil(a,1)))^2);
            nbr_b(xfil(a)+1,yfil(a),1) = bori(xfil(a,1),yfil(a,1));
            nbr_c(xfil(a)+1,yfil(a),1) = cori(xfil(a,1),yfil(a,1));
            nbr_g(xfil(a)+1,yfil(a),1) = gn(xfil(a,1),yfil(a,1));
        end
        %check +x+y
        if yfil(a) ~= W && ind(xfil(a)+1,yfil(a)+1) == 0 && clin(xfil(a)+1,yfil(a)+1) == 0
            pg(xfil(a)+1,yfil(a)+1) = pg(xfil(a)+1,yfil(a)+1) + GP*(sind(cori(xfil(a,1),yfil(a,1))-45)^2+10*cosd(cori(xfil(a,1),yfil(a,1))-45)^2);
            nbr_b(xfil(a)+1,yfil(a)+1,2) = bori(xfil(a,1),yfil(a,1));
            nbr_c(xfil(a)+1,yfil(a)+1,2) = cori(xfil(a,1),yfil(a,1));
            nbr_g(xfil(a)+1,yfil(a)+1,2) = gn(xfil(a,1),yfil(a,1));
        end
        %check +x-y
        if yfil(a) ~= 1 && ind(xfil(a)+1,yfil(a)-1) == 0 && clin(xfil(a)+1,yfil(a)-1) == 0
            pg(xfil(a)+1,yfil(a)-1) = pg(xfil(a)+1,yfil(a)-1) + GP*(sind(cori(xfil(a,1),yfil(a,1))-135)^2+10*cosd(cori(xfil(a,1),yfil(a,1))-135)^2);
            nbr_b(xfil(a)+1,yfil(a)-1,8) = bori(xfil(a,1),yfil(a,1));
            nbr_c(xfil(a)+1,yfil(a)-1,8) = cori(xfil(a,1),yfil(a,1));
            nbr_g(xfil(a)+1,yfil(a)-1,8) = gn(xfil(a,1),yfil(a,1));
        end
    end
    if xfil(a) ~= 1
        %check -x
        if ind(xfil(a)-1,yfil(a)) == 0 && clin(xfil(a)-1,yfil(a)) == 0
            pg(xfil(a)-1,yfil(a)) = pg(xfil(a)-1,yfil(a)) + GP*(10*sind(cori(xfil(a,1),yfil(a,1)))^2+cosd(cori(xfil(a,1),yfil(a,1)))^2);
            nbr_b(xfil(a)-1,yfil(a),5) = bori(xfil(a,1),yfil(a,1));
            nbr_c(xfil(a)-1,yfil(a),5) = cori(xfil(a,1),yfil(a,1));
            nbr_g(xfil(a)-1,yfil(a),5) = gn(xfil(a,1),yfil(a,1));
        end
        %check -x+y
        if yfil(a) ~= W && ind(xfil(a)-1,yfil(a)+1) == 0 && clin(xfil(a)-1,yfil(a)+1) == 0
            pg(xfil(a)-1,yfil(a)+1) = pg(xfil(a)-1,yfil(a)+1) + GP*(sind(cori(xfil(a,1),yfil(a,1))-135)^2+10*cosd(cori(xfil(a,1),yfil(a,1))-135)^2);
            nbr_b(xfil(a)-1,yfil(a)+1,4) = bori(xfil(a,1),yfil(a,1));
            nbr_c(xfil(a)-1,yfil(a)+1,4) = cori(xfil(a,1),yfil(a,1));
            nbr_g(xfil(a)-1,yfil(a)+1,4) = gn(xfil(a,1),yfil(a,1));
        end
        %check -x-y
        if yfil(a) ~= 1 && ind(xfil(a)-1,yfil(a)-1) == 0 && clin(xfil(a)-1,yfil(a)-1) == 0
            pg(xfil(a)-1,yfil(a)-1) = pg(xfil(a)-1,yfil(a)-1) + GP*(sind(cori(xfil(a,1),yfil(a,1))-45)^2+10*cosd(cori(xfil(a,1),yfil(a,1))-45)^2);
            nbr_b(xfil(a)-1,yfil(a)-1,6) = bori(xfil(a,1),yfil(a,1));
            nbr_c(xfil(a)-1,yfil(a)-1,6) = cori(xfil(a,1),yfil(a,1));
            nbr_g(xfil(a)-1,yfil(a)-1,6) = gn(xfil(a,1),yfil(a,1));
        end
    end
    %no need to check +/-y, + is handled by nucleation and there is either
    %only 1 point at the edge (thus no -y), or the line is full and -y will
    %always be filled already.
end

%make struct matrix for holding movies, add first frame
if nameC ~= 0
    ani_pic(t+51) = struct('cdata',[],'colormap',[]);
    ani_pic(1) = im2frame(cori2pic(cori));
end
if nameA ~= 0
    ani_tb(t+51) = struct('cdata',[],'colormap',[]);
    ani_tb(1) = im2frame(GrainBoundMap(bori,cori,gn,1));
end
%% Simulation
for d = 1:t
    if d > 1
        for a = 1:nnew
            if filn(a,1) ~= L
                %check +x
                if ind(filn(a,1)+1,filn(a,2)) == 0 && clin(filn(a,1)+1,filn(a,2)) == 0
                    pg(filn(a,1)+1,filn(a,2)) = pg(filn(a,1)+1,filn(a,2)) + GP*(10*sind(cori(filn(a,1),filn(a,2)))^2+cosd(cori(filn(a,1),filn(a,2)))^2);
                    nbr_b(filn(a,1)+1,filn(a,2),1) = bori(filn(a,1),filn(a,2));
                    nbr_c(filn(a,1)+1,filn(a,2),1) = cori(filn(a,1),filn(a,2));
                    nbr_g(filn(a,1)+1,filn(a,2),1) = gn(filn(a,1),filn(a,2));
                end
                %check +x+y
                if filn(a,2) ~= W && ind(filn(a,1)+1,filn(a,2)+1) == 0 && clin(filn(a,1)+1,filn(a,2)+1) == 0
                    pg(filn(a,1)+1,filn(a,2)+1) = pg(filn(a,1)+1,filn(a,2)+1) + GP*(sind(cori(filn(a,1),filn(a,2))-45)^2+10*cosd(cori(filn(a,1),filn(a,2))-45)^2);
                    nbr_b(filn(a,1)+1,filn(a,2)+1,2) = bori(filn(a,1),filn(a,2));
                    nbr_c(filn(a,1)+1,filn(a,2)+1,2) = cori(filn(a,1),filn(a,2));
                    nbr_g(filn(a,1)+1,filn(a,2)+1,2) = gn(filn(a,1),filn(a,2));
                end
                %check +x-y
                if filn(a,2) ~= 1 && ind(filn(a,1)+1,filn(a,2)-1) == 0 && clin(filn(a,1)+1,filn(a,2)-1) == 0
                    pg(filn(a,1)+1,filn(a,2)-1) = pg(filn(a,1)+1,filn(a,2)-1) + GP*(sind(cori(filn(a,1),filn(a,2))-135)^2+10*cosd(cori(filn(a,1),filn(a,2))-135)^2);
                    nbr_b(filn(a,1)+1,filn(a,2)-1,8) = bori(filn(a,1),filn(a,2));
                    nbr_c(filn(a,1)+1,filn(a,2)-1,8) = cori(filn(a,1),filn(a,2));
                    nbr_g(filn(a,1)+1,filn(a,2)-1,8) = gn(filn(a,1),filn(a,2));
                end
            end
            if filn(a,1) ~= 1
                %check -x
                if ind(filn(a,1)-1,filn(a,2)) == 0 && clin(filn(a,1)-1,filn(a,2)) == 0
                    pg(filn(a,1)-1,filn(a,2)) = pg(filn(a,1)-1,filn(a,2)) + GP*(10*sind(cori(filn(a,1),filn(a,2)))^2+cosd(cori(filn(a,1),filn(a,2)))^2);
                    nbr_b(filn(a,1)-1,filn(a,2),5) = bori(filn(a,1),filn(a,2));
                    nbr_c(filn(a,1)-1,filn(a,2),5) = cori(filn(a,1),filn(a,2));
                    nbr_g(filn(a,1)-1,filn(a,2),5) = gn(filn(a,1),filn(a,2));
                end
                %check -x+y
                if filn(a,2) ~= W && ind(filn(a,1)-1,filn(a,2)+1) == 0 && clin(filn(a,1)-1,filn(a,2)+1) == 0
                    pg(filn(a,1)-1,filn(a,2)+1) = pg(filn(a,1)-1,filn(a,2)+1) + GP*(sind(cori(filn(a,1),filn(a,2))-135)^2+10*cosd(cori(filn(a,1),filn(a,2))-135)^2);
                    nbr_b(filn(a,1)-1,filn(a,2)+1,4) = bori(filn(a,1),filn(a,2));
                    nbr_c(filn(a,1)-1,filn(a,2)+1,4) = cori(filn(a,1),filn(a,2));
                    nbr_g(filn(a,1)-1,filn(a,2)+1,4) = gn(filn(a,1),filn(a,2));
                end
                %check -x-y
                if filn(a,2) ~= 1 && ind(filn(a,1)-1,filn(a,2)-1) == 0 && clin(filn(a,1)-1,filn(a,2)-1) == 0
                    pg(filn(a,1)-1,filn(a,2)-1) = pg(filn(a,1)-1,filn(a,2)-1) + GP*(sind(cori(filn(a,1),filn(a,2))-45)^2+10*cosd(cori(filn(a,1),filn(a,2))-45)^2);
                    nbr_b(filn(a,1)-1,filn(a,2)-1,6) = bori(filn(a,1),filn(a,2));
                    nbr_c(filn(a,1)-1,filn(a,2)-1,6) = cori(filn(a,1),filn(a,2));
                    nbr_g(filn(a,1)-1,filn(a,2)-1,6) = gn(filn(a,1),filn(a,2));
                end
            end
            %check +y
            if filn(a,2) ~= W && ind(filn(a,1),filn(a,2)+1) == 0 && clin(filn(a,1),filn(a,2)+1) == 0
                pg(filn(a,1),filn(a,2)+1) = pg(filn(a,1),filn(a,2)+1) + GP*(sind(cori(filn(a,1),filn(a,2)))^2+10*cosd(cori(filn(a,1),filn(a,2)))^2);
                nbr_b(filn(a,1),filn(a,2)+1,3) = bori(filn(a,1),filn(a,2));
                nbr_c(filn(a,1),filn(a,2)+1,3) = cori(filn(a,1),filn(a,2));
                nbr_g(filn(a,1),filn(a,2)+1,3) = gn(filn(a,1),filn(a,2));
            end
            %check -y
            if filn(a,2) ~= 1 && ind(filn(a,1),filn(a,2)-1) == 0 && clin(filn(a,1),filn(a,2)-1) == 0
                pg(filn(a,1),filn(a,2)-1) = pg(filn(a,1),filn(a,2)-1) + GP*(sind(cori(filn(a,1),filn(a,2)))^2+10*cosd(cori(filn(a,1),filn(a,2)))^2);
                nbr_b(filn(a,1),filn(a,2)-1,7) = bori(filn(a,1),filn(a,2));
                nbr_c(filn(a,1),filn(a,2)-1,7) = cori(filn(a,1),filn(a,2));
                nbr_g(filn(a,1),filn(a,2)-1,7) = gn(filn(a,1),filn(a,2));
            end
        end
    end

    %list coordinates of potential growth/nucleation, find number of them.
    coords = [];
    [coords(:,1),coords(:,2)] = find(pg);
    nc = size(coords,1);

    %make dummy matrices for axis orientations so that pixels added in the
    %same loop don't see each other.
    corin = cori;
    borin = bori;

    %loop through each position, first determine if there is a nucleation,
    %then if there is growth.
    for f = 1:nc
        %check if there is nucleation, if there is set the orientations
        %randomly and move on.
        q1 = rand;
        if q1 <= NP && ismember(coords(f,1),xlin)
            corin(coords(f,1),coords(f,2)) = randi(180);
            borin(coords(f,1),coords(f,2)) = randi(180);
            gc = gc+1;
            gn(coords(f,1),coords(f,2)) = gc;
        else
            %now check if there is growth
            q2 = rand;
            if q2 <= pg(coords(f,1),coords(f,2))
                %if there is, check which neighbors the pixel has and
                %choose one of their orientations to take. Right now
                %which neighbor is selected is purely random.
                nbrt_b = squeeze(nbr_b(coords(f,1),coords(f,2),:));
                nbrt_c = squeeze(nbr_c(coords(f,1),coords(f,2),:));
                nbrt_g = squeeze(nbr_g(coords(f,1),coords(f,2),:));
                nbr_p = (isfinite(nbrt_c)/sum(isfinite(nbrt_c)))';
                q3 = randsample(8,1,true,nbr_p);
                %set the c-axis orientation
                corin(coords(f,1),coords(f,2)) = nbrt_c(q3,1);

                %check for twinning, taking the c-axis orientation into
                %account.
                coeff = TwPcoeff(q3,nbrt_c(q3,1));
                TwP_temp = coeff*(TP);
                q4 = rand;
                if q4 <= TwP_temp
                    %coin toss to see if twinning is +116° or -116°
                    q5 = randi(2);
                    if q5 == 1
                        gc = gc+1;
                        gn(coords(f,1),coords(f,2)) = gc;
                        borin(coords(f,1),coords(f,2)) = nbrt_b(q3,1) + 116;
                        if borin(coords(f,1),coords(f,2)) >= 180
                            borin(coords(f,1),coords(f,2)) = borin(coords(f,1),coords(f,2)) - 180;
                        end
                    else
                        gc = gc+1;
                        gn(coords(f,1),coords(f,2)) = gc;
                        borin(coords(f,1),coords(f,2)) = nbrt_b(q3,1) - 116;
                        if borin(coords(f,1),coords(f,2)) < 0
                            borin(coords(f,1),coords(f,2)) = borin(coords(f,1),coords(f,2)) + 180;
                        end
                    end
                else
                    gn(coords(f,1),coords(f,2)) = nbrt_g(q3,1);
                    borin(coords(f,1),coords(f,2)) = nbrt_b(q3,1);
                end
            end
        end
    end

    %update temporary index
    indn = isfinite(corin);

    %find the number of added pixels this loop and their coordinates by
    %comparing the temporary orientation matrix to the permanent one.
    fil = [];
    filn = [];
    [fil(:,1),fil(:,2)] = find(ind);
    [filn(:,1),filn(:,2)] = find(indn);
    rmv = ismember(filn,fil,'rows');
    %filn is the new coordinates from this loop.
    filn(rmv == 1,:) = [];
    nnew = size(filn,1);

    %remove growth probability from filled coordinates.
    for u = 1:nnew
        pg(filn(u,1),filn(u,2)) = 0;
    end

    %set permanent orientation matrices to new values.
    cori = corin;
    bori = borin;
    ind = indn;

    %record movie frame
    if nameC ~= 0
        ani_pic(d+1) = im2frame(cori2pic(cori));
    end
    if nameA ~= 0
        ani_tb(d+1) = im2frame(GrainBoundMap(bori,cori,gn,1));
    end
end

%fill in the lines that prevent growth purely through random selection of
%neighboring orientations, i.e. no nucleation or twinning, to simulate
%spatial resolution of EBSD where we can't see the lines between individual
%tablets usually
coords = [];
[coords(:,1),coords(:,2)] = find(~ind);
nc = size(coords,1);
corin = cori;
borin = bori;
for m = 1:nc
    %find the orientations of the pixel's neighbors
    nbc = [];
    nbb = [];
    [nbc,nbb] = nbr_ori2(cori,bori,coords(m,1),coords(m,2));
    nbg = nbr_gn(gn,coords(m,1),coords(m,2));
    nbc = rmmissing(nbc);
    nbb = rmmissing(nbb);
    nbg = rmmissing(nbg);

    %randomly choose one, set the value in the temp matrix
    if isempty(nbc) == 0
    qn = randi(size(nbc,2));
        corin(coords(m,1),coords(m,2)) = nbc(qn);
        borin(coords(m,1),coords(m,2)) = nbb(qn);
        gn(coords(m,1),coords(m,2)) = nbg(qn);
    end
end
    
%set the temp orientations into the final matrices
bori = borin;
cori = corin;

%record movie frame for the filled in 
if nameC ~= 0
    ani_pic(d+1:d+51) = im2frame(cori2pic(cori));
end
if nameA ~= 0
    ani_tb(d+1:d+51) = im2frame(GrainBoundMap(bori,cori,gn,1));
end
%% Twin boundary map
%map the twin boundaries in the center slice of the simulation. Currently
%assigns each pixel the highest twin multiple it has adjacent to it
%comparing in the positive direction.
tim = GrainBoundMap(bori,cori,gn,1);
imshow(tim);
figure()
%% Generate c-axis PIC map
pim = cori2pic(cori);
imshow(pim)
%% Make .mp4 movie files
if nameC ~= 0
    v1 = VideoWriter(nameC,'MPEG-4');
    open(v1)
    writeVideo(v1,ani_pic);
    close(v1)
end

if nameA ~= 0
    v2 = VideoWriter(nameA,'MPEG-4');
    open(v2)
    writeVideo(v2,ani_tb);
    close(v2)
end