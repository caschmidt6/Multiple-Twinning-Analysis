% crystal symmetry
% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('mmm', [5.7 5 8], 'mineral', 'Aragonite', 'color', [0.53 0.81 0.98]),...
  'notIndexed'};

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

%% Specify File Names

% path to files
pname = 'C:\Users\Connor Schmidt\Desktop\Raw EBSD data\Hrf3-11 nacre\Hrf3_11area8.ctf';

% which files to be imported CHANGE THIS NAME TO CHANGE FILES
fname = [pname '\Hrf3_11area8.ctf'];

%% Import the Data

% create an EBSD variable containing the data
eHrf8 = EBSD.load(fname,CS,'interface','ctf',...
  'convertSpatial2EulerReferenceFrame');

%Find the grains, eliminating very small ones
[gGeo3,eHrf8.grainId,eHrf8.mis2mean] = calcGrains(eHrf8,'alpha',2,'angle',5*degree);
eHrf8(gGeo3(gGeo3.grainSize<5)) = [];
[gGeo3,eHrf8.grainId,eHrf8.mis2mean] = calcGrains(eHrf8,'alpha',2,'angle',5*degree);