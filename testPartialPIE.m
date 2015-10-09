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
del = 4E-8;
Lc = 10E-6;
TotalSteps = 30;
suffix='exp';

%Pos = dlmread(fullfile(src,'Position.txt'));
%load(fullfile(src,'PtychoData.mat'));
if exist('Step','var') 
    if Step>1E-5
        Step = Step*1E-6;
    end
else
    Step = 1E-6;
end
if exist('Lc','var')
    if Lc>1
        Lc = Lc*1E-6;
    end
else
    Lc = Inf;
end
PWD='/home/bchen/ASync019_scratch/BoChen';
fprintf(1,'***************************************\n');
fprintf(1,'** Reconstruction Ptychography **\n');
fprintf(1,'***************************************\n');
fprintf(1,'Parameters\n');
fprintf(1,'del = %E\n',del);
fprintf(1,'E = %d\n',E);
fprintf(1,'z1 = %E\n',z1);
fprintf(1,'Lc= %E\n',Lc);
fprintf(1,'Step= %E\n',Step);
fprintf(1,'***************************************\n');

% data_src = sprintf('CDI/Scan%g/Lc%d/%s',Step*1E6,...
%     round(Lc*1E6),suffix);
% fullname = fullfile(PWD,data_src,'SimData.h5');
fullname = fullfile('/home/bchen/ASync019_scratch/BoChen/CDI/Scan2.5/Lc8/ideal/SimData.h5');
Ie = hdf5read(fullname,'/data/intensity');
Pos = hdf5read(fullname,'/data/position');
Probe = hdf5read(fullname,'/data/probe_real')+...
    1j*hdf5read(fullname,'/data/probe_imag');
M = size(Probe,1);
xpos = -Pos(:,1);
ypos = -Pos(:,2);
if noiseLevel<0 
    recon = 'ePIE_recon_partial';
else
    recon=sprintf('ePIE_recon_Noise_%gpct_partial',noiseLevel);
    Ie = Ie + noiseLevel/100*mean(Ie(:))*randn(size(Ie));
    Ie(Ie<0) = 0;
end
if AID>=0
    recon = sprintf('%s/%03d',recon,AID);
end
dest = fullfile(data_src,recon);
[V,D]=CalcModes(Lc,del,M,numModes,numModes);
partialePIE;
