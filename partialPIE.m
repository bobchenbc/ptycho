function results = partialPIE(params)
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
%
params = setDefault(params);

Power = sum(sum(params.Ie,1),2);
Power = squeeze(Power);

[M,N,numpts] = size(params.Ie);
numModes = params.numModes;
xpos = params.xpos;
ypos = params.ypos;
obs = params.obs;


mag = sqrt(params.Ie);
for k=1:numpts
    mag(:,:,k) = fftshift(mag(:,:,k));
end

probe = zeros(M,N,numModes);

params.Probe = params.Probe./sqrt(sum(abs(params.Probe(:).^2)));
params.Probe = fftshift(ifft2(fftshift(params.Probe)));
for k = 1:numModes
    probe(:,:,k) = fftshift(fft2(fftshift(params.Probe.*params.V(:,:,k))));
    s = sqrt(sum(sum(abs(probe(:,:,k)).^2)));
    if s>0
        probe = probe/s;
    end
end
probe = probe*(sum(abs(params.Probe(:)).^2)/sum(abs(probe(:)).^2));
params=rmfield(params,{'xpos','ypos','obs','Probe','Ie'});


%probe = repmat(Probe,[1,1,numModes])./numModes;
psix = complex(zeros(M,N,numModes),0);
po = complex(zeros(M,N,numModes),0);

err = 1E3;
Err = zeros(params.TotalSteps,1);
step = 1;

if params.GPU
    probe = gdouble(probe);
    psix = gdouble(psix);
    po = gdouble(po);
    obs = gdouble(obs);
    Power = gdouble(Power);
end

while step<=params.TotalSteps
    if params.GPU
        err = gdouble(0);
    end
    for j=1:numpts
        indy = (1:M)+ypos(j);
        indx = (1:N)+xpos(j);
        % for can be replaced with parfor
        for k=1:numModes
            po(:,:,k) = probe(:,:,k).*obs(indy,indx);
            psix(:,:,k) = fft2(po(:,:,k));
        end
        % calculated magnitude
        magc = sqrt(sum(abs(psix).^2,3));
        err = err + sum(sum((magc-mag(:,:,j)).^2))./Power(j);

        for k=1:numModes
            tmp = psix(:,:,k)./magc;
            tmp(isnan(tmp)|isinf(tmp))=0;
            psix(:,:,k) = mag(:,:,j).*tmp;
            psix(:,:,k) = ifft2(psix(:,:,k));
            df = psix(:,:,k) - po(:,:,k);
            mx = max(max(abs(probe(:,:,k))));
            obs(indy,indx) = obs(indy,indx) + conj(probe(:,:,k))./mx.^2.*df;
            if step>params.updateProbeSteps
                mx = max(max(abs(obs(indy,indx))));
                probe(:,:,k) = probe(:,:,k) + conj(obs(indy,indx))/mx.^2.*df;
            end
        end
    end

    %obs=custom_constraint(obs,'forceunity');
    if params.showFigure
        if isa(obs,'gpuArray')
            showmat(angle(gather(obs)));
        else
            showmat(angle(obs));
        end
        drawnow;
    end
    err = err/numpts;


    if isa(obs,'gpuArray')
        this_err = gather(err);
    else
        this_err=err;
    end
    Err(step) = this_err;
    fprintf(1,'step=%03d, err=%.3f\n',step,this_err);

    if mod(step, params.outputSteps) ==0
        if isa(obs,'gpuArray')
            im_mag = gather(abs(obs));
            im_phs = gather(angle(obs));
        else
            im_mag = abs(obs);
            im_phs = angle(obs);
        end

        im = mat2img(im_mag,1);
        imwrite(im,fullfile(params.dest,sprintf('mag_step_%04d.png',step)));
        im = mat2img(im_phs,1);
        imwrite(im,fullfile(params.dest,sprintf('phs_step_%04d.png',step)));

    end
    step = step+1;

end

results.obs = gather(obs);
results.Err = gather(Err);
results.probe = gather(probe);
%results.pos = [xpos,ypos];
writeData(results,params);
end

function writeData(results,params)
    if strcmpi(params.saveType,'mat')
        fullname = fullfile(params.dest,'result.mat');
        fprintf(1,'Writing results to %s...',fullname);
        save(fullname,'results');
        fprintf(1,'Done!\n');
    end
    if strcmpi(params.saveType,'h5')
        fullname = fullfile(params.dest,'result.h5');

        fprintf(1,'Writing /result/obs_real...');
        hdf5write(fullname,'/result/obs_real',real(results.obs));
        fprintf(1,'Done!\nWriting /result/obs_imag...');
        hdf5write(fullname,'/result/obs_imag',imag(results.obs),'WriteMode','append');
        fprintf(1,'Done!\nWriting /result/probe_real...');
        hdf5write(fullname,'/result/probe_real',real(results.probe),'WriteMode','append');
        fprintf(1,'Done!\nWriting /result/probe_imag...');
        hdf5write(fullname,'/result/probe_imag',imag(results.probe),'WriteMode','append');
        fprintf(1,'Done!\nWriting /result/error_metric...');
        hdf5write(fullname,'/result/error_metric',results.Err,'WriteMode','append');
        fprintf(1,'Done!\n');
    end
end

function params=setDefault(params)

    % number of spatially coherent modes
    if ~isfield(params,'numModes')
        params.numModes = 1;
    end

    % outputSteps: steps to output temporal results, default 20
    if ~isfield(params,'outputSteps')
        params.outputSteps = 20;
    end

    % total iteration steps
    if ~isfield(params,'TotalSteps')
        params.TotalSteps = 50;
    end

    if ~isfield(params,'method')
        params.method = 'ePIE';
    end

    % number of extra pixels in obs, it is useful in position correction
    if ~isfield(params,'relax')
        params.relax = 0;
    end



    if ~isfield(params,'Ie')
        error('Please provide experimental diffraaction intensity(Ie)');
    end
    if ~isfield(params,'xpos')||~isfield(params,'ypos')
        error('Please provide the scanning positions in pixel number(xpos,ypos)');
    end
    if ~isfield(params,'Probe')
        error('Please provide the initial probe function(Probe)');
    end


    % pixel size in the sample plane
    if ~isfield(params,'del')
        error('Please provide pixel size in the sample plane(del)');
    end
    [M,N,~] = size(params.Ie);
    %
    % initial guess of the coherence length
    if ~isfield(params,'Lc')
        params.Lc = 10E-6;
    elseif isinf(params.Lc)
        params.Lc = M/2*params.del;
    end


    [V,D]=CalcModes(params.Lc,params.del,M,params.numModes);

    for k =1:length(D)
        V(:,:,k) = V(:,:,k).*D(k);
    end
    params.V = V;

    % updateProbeSteps: steps when probe begin to update
    if ~isfield(params,'updateProbeSteps')
        params.updateProbeSteps = 10;
    end

    if ~isfield(params,'showFigure')
        params.showFigure=false;
    end

    % dest: Destination folder to write results, default 'result[xxx]'
    if ~isfield(params,'dest')
        n = 1;
        params.dest = sprintf('result%03d',n);
        while exist(params.dest,'dir')
            n = n+1;
            params.dest = sprintf('result%03d',n);
        end
    end
    if ~exist(params.dest,'dir')
        mkdir(params.dest);
    end

    if ~isfield(params,'initType')
        params.initType= 'flat';
    end

    if ~isfield(params,'saveType')
        params.saveType='mat';
    end

    if ~isfield(params,'obs')
        params.xpos = params.xpos - min(params.xpos);
        params.ypos = params.ypos - min(params.ypos);
        params.xpos = params.xpos + round(params.relax/2);
        params.ypos = params.ypos + round(params.relax/2);

        ex = round(max(params.xpos)+params.relax/2);
        ey = round(max(params.ypos)+params.relax/2);
        switch params.initType
            case 'flat'
                params.obs = complex(ones(M+ey,N+ex),0);
            otherwise
                params.obs = complex(rand(M+ey,N+ex),0);
        end
    end

    if ~isfield(params,'GPU')
        params.GPU=false;
    end
end
