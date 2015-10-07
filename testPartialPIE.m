numModes = 5;
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

Pos = dlmread(fullfile(src,'Position.txt'));
load(fullfile(src,'PtychoData.mat'));

D = 10E-6; % diameter of probe
R = D/2;
del = 4E-8;

NF = R.^2/(lambda*z1);

[xi,yi] = meshgrid((1:M)-round(M/2+0.5));
z2 = (xi).^2+(yi).^2;
Probe = double(z2<=(R/del).^2);
% u = fftshift(fft2(fftshift(Probe)));
% fac = 1i*pi*(f+fd)/(f*fd)*(df^2/lambda);
% phase = exp(fac*(xi.^2+yi.^2));
% u = fftshift(fft2(fftshift(phase.*u)));
% fac = -1i*pi*fs/(fd*(fd-fs))*(ps^2/lambda);
% phase = exp(fac*(xi.^2+yi.^2));
% Probe = fftshift(ifft2(fftshift(u.*phase)));   
if NF>1
    Probe = prop_free_ff(Probe,lambda,z1,del);
else
    Probe = prop_free_nf(Probe,lambda,z1,del);
end


partialePIE;

