numModes = 3;
outputSteps = 20;
updateProbeSteps = 5;
relax=0;
threshold = 1E-5;
dest = 'resulttest';
src = './';
z1 = 1E-3;
E = 1400;
lambda = 12.4E-7/E;
M = 1024;

%Pos = dlmread(fullfile(src,'Position.txt'));
%load(fullfile(src,'PtychoData.mat'));

data_src = '/home/bchen/ASync019_scratch/BoChen/CDI/Scan2/Lc6/ideal/';
fullname = fullfile(data_src,'SimData.h5');
Ie = hdf5read(fullname,'/data/intensity');
Pos = hdf5read(fullname,'/data/position');
Probe = hdf5read(fullname,'/data/probe_real')+...
    1j*hdf5read(fullname,'/data/probe_imag');

xpos = -Pos(:,1);
ypos = -Pos(:,2);
partialePIE;
