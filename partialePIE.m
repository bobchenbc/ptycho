% multi-modes ptychography reconstruction algorithm
% parameters needed
% numModes: number of modes, default 1
% outputSteps: steps to output temporal results, default 20
% relax: number of extra pixels to add in obs, default 0
% threshold: the minimum value of the calculated magnitude. Any value less than
% threhold will be set to zero, default value = 1E-5;
% updateProbeSteps: steps when probe begin to update
% less than threshold is equal to threhold
% dest: Destination folder to write results, default 'result[xxx]'
%
% Probe
% Pos
% Ie
% 



% number of spatially coherent modes
if ~exist('numModes','var')
    numModes = 1;
end

% outputSteps: steps to output temporal results, default 20
if ~exist('outputSteps','var')
    outputSteps = 20;
end

% number of extra pixels in obs, it is useful in position correction
if ~exist('relax','var')
    relax = 0;    
end

% threshold: the minimum value of the calculated magnitude. Any value 
if ~exist('threshold','var')
    threshold = 1E-5;
end

% updateProbeSteps: steps when probe begin to update
if ~exist('updateProbeSteps','var')
    updateProbeSteps = 10;
end

% dest: Destination folder to write results, default 'result[xxx]' 
if ~exist('dest','var')
    n = 1;
    dest = sprintf('result%03d',n);
    while exist(dest,'dir')
        n = n+1;
         dest = sprintf('result%03d',n);
    end
end
if ~exist(dest,'dir')
    mkdir(dest);
end

if ~exist('Ie','var')
    error('Please provide experimental diffraaction intensity(Ie)');
end

if ~exist('Pos','var')
    error('Please provide the scanning positions in pixel number(Pos)');
end

if ~exist('Probe','var')
    error('Please provide the initial probe function(Probe)');
end


Power = sum(sum(Ie,1),2);
Power = squeeze(Power);
maxPower = max(Power);
[M,N,numpts] = size(Ie);

for k=1:numpts
    Ie(:,:,k) = fftshift(Ie(:,:,k)); %#ok<SAGROW>
end
mag = sqrt(Ie);

Probe = Probe.*sqrt(maxPower/(M*N*sum(abs(Power(:).^2))));
probe = repmat(Probe,[1,1,numModes])./sqrt(numModes);
psix = complex(zeros(M,N,numModes),0);
po = complex(zeros(M,N,numModes),0);

xpos = Pos(:,1);
ypos = Pos(:,2);
xpos = xpos - min(xpos);
ypos = ypos - min(ypos);


xpos = xpos + round(relax/2);
ypos = ypos + round(relax/2);
ex = round(max(xpos)+relax/2);
ey = round(max(ypos)+relax/2);

if ~exist('obs','var')
    obs = complex(rand(M+ey,N+ey),0);
end


omega = 1E3;
step = 1;

while omega > 1E-2
    err = 0;
    for k=1:numModes
        for j=1:numpts
            indy = (1:M)+ypos(j);
            indx = (1:N)+xpos(j);
            po(:,:,k) = probe(:,:,k).*obs(indy,indx);
            psix(:,:,k) = fft2(po(:,:,k));
        end
        % calculated magnitude
        magc = sqrt(sum(abs(psix).^2,3));
        err = err + sum(sum((magc-mag(:,:,j)).^2))./Power(j);
        fprintf(1,'err=%g\n',err);
        magc(magc<threshold) = threshold;
        magc = mag(:,:,k)./magc;
        
        % for can be replaced with parfor
        for k=1:numModes 
            psix(:,:,k) = magc.*psix(:,:,k);%exp(1i*angle(psix(:,:,n)));
            psix(:,:,k) = ifft2(psix(:,:,k));
            df = psix(:,:,k) - po(:,:,k);
            mx = max(max(abs(probe(:,:,k))));
            obs(indy,indx) = obs(indy,indx) + conj(probe(:,:,k))./mx.^2.*df;
        end
        if step>updateProbeSteps
            for k =1:numModes
                df = psix(:,:,k) - probe(:,:,k).*obs(indy,indx);
                mx = max(max(abs(obs(indy,indx))));
                po(:,:,k) = po(:,:,k) + conj(obs(indy,indx))/mx.^2.*df;
            end
        end
                
    end
    
    fprintf(1,'step=%03d, err=%.3f\n',step,err);
    step = step+1;
end
