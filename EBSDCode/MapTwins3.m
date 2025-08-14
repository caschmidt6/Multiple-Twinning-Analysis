function [Twins,hist] = MapTwins3(ebsd,grains,maxT,err)
%"ebsd" should be an EBSD variable, "grains" should be a grains2d variable,
%"maxT" is the highest twin multiple to be labeled, and "err" is the
%allowed range around the true twin misorientation to be considered that
%twin multiple (i.e. err = 1 means twin 1 can have a misorientation
%anywhere from 62.8° to 64.8°).
%"Twins" is the list of twin boundaries (1 line per boundary segment, nan
%for not a twin, 0 for co-oriented in c but not a twin, or numbered with
%the twin multiple). "hist" is the number of each twin multiple in the
%given area.

%identify twin boundary segments from inputs
Twins = FindTwins(grains,maxT,err);

%
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