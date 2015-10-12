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





    Power = sum(sum(Ie,1),2);
    Power = squeeze(Power);
    maxPower = max(Power);
    [M,N,numpts] = size(Ie);

    for k=1:numpts
        Ie(:,:,k) = fftshift(Ie(:,:,k)); %#ok<SAGROW>
    end
    mag = sqrt(Ie);
    magc = ones(size(mag));

    probe = zeros(M,N,numModes);
    mx = max(Probe(:));
    Probe = Probe./sqrt(sum(abs(Probe(:).^2)));
    Probe = fftshift(ifft2(fftshift(Probe)));
    for k = 1:numModes
        probe(:,:,k) = fftshift(fft2(fftshift(Probe.*V(:,:,k)*D(k))));
        probe = probe/sqrt(sum(sum(abs(probe(:,:,k)).^2)));
    end

    probe = probe*(sum(abs(Probe(:)).^2)/sum(abs(probe(:)).^2));




    %probe = repmat(Probe,[1,1,numModes])./numModes;
    psix = complex(zeros(M,N,numModes),0);
    po = complex(zeros(M,N,numModes),0);

    xpos = xpos - min(xpos);
    ypos = ypos - min(ypos);


    xpos = xpos + round(relax/2);
    ypos = ypos + round(relax/2);
    ex = round(max(xpos)+relax/2);
    ey = round(max(ypos)+relax/2);

    if ~isfield(params,'obs')
        obs = complex(ones(M+ey,N+ey),0);
    end


    omega = 1E3;
    step = 1;

    while omega > 1E-2 && step<=TotalSteps
        err = 0;
        for j=1:numpts
            indy = (1:M)+ypos(j);
            indx = (1:N)+xpos(j);
            for k=1:numModes
                po(:,:,k) = probe(:,:,k).*obs(indy,indx);
                psix(:,:,k) = fft2(po(:,:,k));
            end
            % calculated magnitude
            magc(:,:,j) = sqrt(sum(abs(psix).^2,3));
            err = err + sum(sum((magc(:,:,j)-mag(:,:,j)).^2))./Power(j);
            %fprintf(1,'err=%g\n',err);
            %end

            % for can be replaced with parfor
            for k=1:numModes 
                %for j = 1:numpts
                %indy = (1:M)+ypos(j);
                %indx = (1:N)+xpos(j);
                %po(:,:,k) = probe(:,:,k).*obs(indy,indx);
                %psix(:,:,k) = fft2(po(:,:,k));
                psix(:,:,k) = mag(:,:,j).*(psix(:,:,k)./magc(:,:,j));
                %psix(:,:,k) = mag(:,:,j).*exp(1i*angle(psix(:,:,k)));
                psix(:,:,k) = ifft2(psix(:,:,k));
                df = psix(:,:,k) - po(:,:,k);
                mx = max(max(abs(probe(:,:,k))));
                obs(indy,indx) = obs(indy,indx) + conj(probe(:,:,k))./mx.^2.*df;
                if step>updateProbeSteps
                    mx = max(max(abs(obs(indy,indx))));
                    probe(:,:,k) = probe(:,:,k) + conj(obs(indy,indx))/mx.^2.*df;
                end

            end
        end

        %obs=custom_constraint(obs,'forceunity');
        showmat(angle(obs));
        drawnow;
        err = err/numpts;
        fprintf(1,'step=%03d, err=%.3f\n',step,err);
        if mod(step, outputSteps) ==0
            mag = abs(obs);
            im = mat2img(mag,1);
            imwrite(im,fullfile(dest,sprintf('mag_step_%04d.png',step)));
            phs = angle(obs);
            im = mat2img(phs,1);
            imwrite(im,fullfile(dest,sprintf('phs_step_%04d.png',step)));

        end
        step = step+1;
    end
    fullname = fullfile(dest,'result.h5');
    fprintf(1,'Writing /result/obs_real...');
    hdf5write(fullname,'/result/obs_real',real(obs));
    fprintf(1,'Done!\nWriting /result/obs_imag...');
    hdf5write(fullname,'/result/obs_imag',imag(obs),'WriteMode','append');
    fprintf(1,'Done!\nWriting /result/probe_real...');
    hdf5write(fullname,'/result/probe_real',real(probe),'WriteMode','append');
    fprintf(1,'Done!\nWriting /result/probe_imag...');
    hdf5write(fullname,'/result/probe_imag',imag(probe),'WriteMode','append');
    fprintf(1,'Done!\nWriting /result/error_metric...');
    hdf5write(fullname,'/result/error_metric',Err,'WriteMode','append');
    fprintf(1,'Done!\n');

function params=setDefault(params)

    % number of spatially coherent modes
    if ~isfield(params,'numModes')
        params.numModes = 1;
    end

    % outputSteps: steps to output temporal results, default 20
    if ~isfield(params,'outputSteps')
        params.outputSteps = 20;
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

    % threshold: the minimum value of the calculated magnitude. Any value 
    if ~isfield(params,'threshold')
        parmas.threshold = mean(Ie(:))*1E-2;
    end

    % updateProbeSteps: steps when probe begin to update
    if ~isfield(params,'updateProbeSteps')
        params.updateProbeSteps = 10;
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

end
