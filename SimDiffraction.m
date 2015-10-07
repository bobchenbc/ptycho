Lcs = 10E-6;
M = 1024;
pix = 13.5E-6;
del = 4E-8;
E = 1400;
z1 = 1E-3;%1E-3; % distance between sample and pinhole
steps=2E-6;
%d =10E-6; % diameter of the pinhole
D = 10E-6; % diameter of probe

R = D/2;
lambda = 12.4E-7/E;
z = M*del*pix/lambda; %
NFS = R^2/(lambda*z);%(NFS=0.0452, NF for a 10um pinhole)
%NF = (M*del).^2/(2*lambda*z);%(NF=1.517 NF for the whole matrix)

g = readcplx('test_biological_1.bin');
%g = hole([M,M], [], round(d/2/del));
%psf = fspecial('average',4);
%g = imfilter(double(g),psf);


C = round(M/2+0.5);

[xi,yi] = meshgrid((1:M)-C);


map = exp(1i*pi*((xi*del).^2+(yi*del).^2)/(lambda*z));

p.Sample2CCD=z;
p.CCDGain=1;
p.NumFrames = 1;
p.PixelSize = pix;
p.RealData = false;
p.BeamStopArea = 0;


p.Energy = E;

fprintf(1,'***************************************\n');
fprintf(1,'** Simulate Ptychography **\n');
fprintf(1,'***************************************\n');

z2 = (xi).^2+(yi).^2;
probe = double(z2<=(R/del).^2);
%psf = fspecial('gaussian', 5,2);
%probe = imfilter(probe,psf);

NF = R.^2/(lambda*z1);

if NF<1
    probe = prop_free_ff(probe,lambda,z1,del);
else
    probe = prop_free_nf(probe,lambda,z1,del);
end
probe = probe./max(probe(:));
%   probe(probe<0.013)=0;
if p.RealData
    suffix = 'roundoff';
else
    suffix = 'raw';
end

if z1>0
    suffix = [suffix,'_prop'];
else
    suffix = [suffix,'_noprop'];
end

for step = steps
    [ns,ms]= meshgrid((-3:3)*step);
    ms = ms(:);
    ns = ns(:);
    %ms = ms + (4*rand(size(ms))-2)*del;
    %ns = ns + (4*rand(size(ns))-2)*del;
    ms = round(ms/del);
    ns = round(ns/del);

    I = zeros(M,M,length(ms));
    for LK = 1:length(Lcs)
        Lc = Lcs(LK);
        %dirname = sprintf('Scan%gum/D%dLc%d/%s',step*1E6,...
        %    round(D*1E6),round(Lc*1E6),suffix);
        dirname = './';
        if ~exist(dirname,'dir')
            mkdir(dirname);
        end
        p.CoherenceLength =Lc;
        for k = 1:length(ms)        
            x0 = ms(k);
            y0 = ns(k);

            A = probe .* circshift(g,[y0,x0]);
            p.g = A;

            Ie = PartialSim(p);

            subplot(1,2,1);
            showmat(A);
            axis image;
            subplot(1,2,2);
            imagesc(log(Ie));
            axis image;
            drawnow;
            %pause(1);

            I(:,:,k) = Ie;
            %             if k== round(length(ms)/2+0.5)
            %                 p.g = g;
            %                 Ie = PartialSim(p);
            %                 
            %                 fullname = fullfile(dirname,'SingleData.dbin2');
            %                 writedbin2(fullname,Ie);
            %             end
            %             fullname = fullfile(dirname,sprintf('data%02d.dbin',k));
            %             writedbin(fullname,Ie);

        end
        fullname = fullfile(dirname,'PtychoData.mat');
        save(fullname,'I');
        dlmwrite(fullfile(dirname,'Position.txt'),round([ms,ns]));
    end
end
scsz = get(0,'ScreenSize');
if sum(scsz(3:4))<3
    quit force
end
