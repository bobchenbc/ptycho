params.numModes = 3;
params.outputSteps = 20;
params.updateProbeSteps = 5;
params.relax=0;
%params.src = '/home/bchen/ASync019_scratch/BoChen/CDI/Scan2.5/Lc8/ideal/';
params.del = 4E-8;
params.Lc = 10E-6;
params.TotalSteps = 100;
PWD='/home/bchen/ASync019_scratch/BoChen';
suffix='exp';
%noiseLevel=1;
%AID=0;

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
params.Lc = Lc;

if ~exist('noiseLevel','var')
    noiseLevel = 0;
end
if ~exist('AID','var')
    AID= -1;
end
fprintf(1,'***************************************\n');
fprintf(1,'** Reconstruction Ptychography **\n');
fprintf(1,'***************************************\n');
fprintf(1,'Parameters\n');
fprintf(1,'params.del = %E\n',params.del);
fprintf(1,'params.Lc= %E\n',params.Lc);
fprintf(1,'params.Step= %E\n',Step);
fprintf(1,'***************************************\n');

params.src = sprintf('/home/bchen/ASync019_scratch/BoChen/CDI/Scan%g/Lc%d/%s',Step*1E6,...
    round(Lc*1E6),suffix);
fullname = fullfile(params.src,'SimData.h5');
params.Ie = hdf5read(fullname,'/data/intensity');
Pos = hdf5read(fullname,'/data/position');
params.Probe = hdf5read(fullname,'/data/probe_real')+...
    1j*hdf5read(fullname,'/data/probe_imag');
M = size(params.Probe,1);
params.xpos = -Pos(:,1);
params.ypos = -Pos(:,2);
if noiseLevel<0 
    recon = 'ePIE_recon_partial';
else
    recon=sprintf('ePIE_recon_Noise_%gpct_partial',noiseLevel);
    params.Ie = params.Ie + noiseLevel/100*mean(params.Ie(:))*randn(size(params.Ie));
    params.Ie(params.Ie<0) = 0;
end
if AID>=0
    recon = sprintf('%s/%03d',recon,AID);
end
if AID==0
    params.initType='flat';
else
    params.initType='rand';
end
params.dest = fullfile(params.src,recon);
params.GPU=true;

results=partialPIE(params);
if AID==0
    M = 325;
    N = 314;
    t1 = cropmat(results.obs,M,N);
    writecplx2img(fullfile(params.dest,'final.png'),t1,true);
end
