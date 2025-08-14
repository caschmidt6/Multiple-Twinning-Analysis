function [Twins,hist] = MapTwins3(ebsd,grains,maxT,err)
%plot twin boundaries in a color-coded scale from 0x to 6x twins, labelling
%only grain boundaries with c-axis misorientation less than or equal to the
%threshold, and only labelling a twin if the corresponding twin
%multiplicity appears as a peak on the histogram (i.e. the corresponding
%b-axis misorientation is more common than the surrounding values).

Twins = FindTwins(grains,maxT,err);

gbs = grains.boundary.grainId;
[~,uq,~] = unique(gbs,'rows');
hist = rmmissing(Twins(uq));

plot(ebsd,ebsd.bc)
hold on
plot(grains.boundary(Twins == 0),'linewidth',2)
plot(grains.boundary(Twins == 1),'linewidth',2,'linecolor','red')
if maxT >= 2
    plot(grains.boundary(Twins == 2),'linewidth',2,'linecolor',[1 0.5 0])
end
if maxT >= 3
    plot(grains.boundary(Twins == 3),'linewidth',2,'linecolor',[1 1 0])
end
if maxT >= 4
    plot(grains.boundary(Twins == 4),'linewidth',2,'linecolor',[0 1 0])
end
if maxT >= 5
    plot(grains.boundary(Twins == 5),'linewidth',2,'linecolor',[0 1 0.71])
end
if maxT >= 6
    plot(grains.boundary(Twins == 6),'linewidth',2,'linecolor',[0 1 1])
end
if maxT >= 7
    plot(grains.boundary(Twins == 7),'linewidth',2,'linecolor',[0 0.5 1])
end
if maxT >= 8
    plot(grains.boundary(Twins == 8),'linewidth',2,'linecolor',[0 0 1])
end
if maxT >= 9
    plot(grains.boundary(Twins == 9),'linewidth',2,'linecolor',[0.5 0 1])
end