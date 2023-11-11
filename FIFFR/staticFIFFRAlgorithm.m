%% FINAL FIFFR APPROACH
%  PROGRESSIVE-RELAXED FILTERING
%  cnc Mar 21, 2023

%% LOAD DATA
data = load('SparseSampledButterflyData.mat'); 
Nx = data.Nx; % x image dimensions
Ny = data.Ny; % y image dimensions
Original = data.Original; % fully sampled image
ffarray = data.ffarray; % array of the fill factors assigned with each image
fillLIST = data.fillLIST; 
staticFIG = data.staticFIG; % array of sparsely sampled versions of the original image, with fill factors corresponding to ffarray
staticIDX = data.staticIDX; % indices corresponding to each image in staticFIG, sampled pixels = 1, unsampled pixels = 0
%% Parameters

% units
nm = 1e-9; 
um = 1e-6; 
mm = 1e-3; 
cm = 1e-2; 

% pixel size in [m]
pixelsize = 2.49e-6; 

% colormap
vec = [0.25; 0];
hex = ['#ffffff'; '#000000'];
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
N = 128;
%N = size(get(gcf,'colormap'),1) % size of the current colormap
map = interp1(vec,raw,linspace(0, 100,N),'pchip');
map = map./(max(map));

% functions
normalize2D =@(M) (M-min(min(M)))/(max(max(M))-min(min(M))); % normalize an image
normalizeM = @(M) (M-min(M))/(max(M)-min(M)); % normalize a 2D vector
sigmoid =@(a, x, c, maxS, startS) (maxS-startS).*normalizeM(1./(1 + exp(-a.*(x-c))))+startS;  % sigmoid function


% selections
% select which image you would like to cycle through FIFFR
imind = 8; % IMage INDex; 1 is the least sparsely sampled image, 10 is the most sparse image
cyclenumber = 150; % number of iterations to cycle FIFFR
replvalue = 0; % Value used to intialize unsampled pixels

%%%%%%%%%%%%%%%
%% Run FIFFR %%
%%%%%%%%%%%%%%%

% initializations
Irecon = staticFIG{imind}; % grab the image
Irecon = normalize2D(Irecon); % normalize
index = (staticFIG{imind}>0); % grab the indices associated with sampled pixels
NOriginal = normalize2D(Original); % grab the original image, fully sampled

% replace all null/NaN values in the sparsely sampled image with the
% replacement value
Irecon(isnan(Irecon)) = replvalue; 

% Variables to track: normalized
FirstImage = Irecon;
VeryFirstImage = Irecon; 

% mask initial params
[Ny, Nx] = size(Irecon); % image size
xgrid = linspace(0, 1e-3, Nx); % x vector
dx = xgrid(2)-xgrid(1); % x sample spacing
Dx = max(xgrid)-min(xgrid); 
dfx = 1/Dx; 
fx = dfx*([1:Nx]-Nx/2-1); % fx vector
ygrid = linspace(0, max(xgrid)*(Ny/Nx), Ny); % y vector
dy = ygrid(2)-ygrid(1); 
Dy = max(ygrid)-min(ygrid); 
dfy = 1/Dy; 
fy = dfy*([1:Ny]-Ny/2-1); % fy vector
[Xt, Yt] = meshgrid(xgrid,ygrid); 
Rt = sqrt(Xt.^2+Yt.^2); 
[FXt, FYt] = meshgrid(fx, fy); 
FRt = sqrt(FXt.^2 + FYt.^2); 

NA = 0.13; % numerical aperture
Mag = 0.5; % f2/f1 (magnification)
ftube = 15*mm; % this is f2
n = 1; % index of refraction
lambda = 785*nm; % wavelength
fobj = ftube/Mag; % this is f1

firstMax = max(fx)/(NA/(n*lambda)); % maximum spatial frequency in x dimension
secondMax = max(fy)/(NA/(n*lambda)); % maximum spatial frequency in y dimension
maxkmask = 2*max(firstMax, secondMax);  % maximum spatial frequency 

% parameters for the thresholding function
a = 0.3; 
sigshift = 9.5;  
maxthresh = maxkmask;
minthresh = 0.1;

% create threshold function
sx = linspace(-10, 10,cyclenumber); 
pthresh = sigmoid(a, sx, sigshift, maxthresh, minthresh); 

tic
for i = 1:length(pthresh)
    FTimage = fftshift(fft2(fftshift(Irecon))); % FFT of sparse image for i=1; FT of estimate for i>1
    Filt = FTimage; % store FT 
    [mask] = generateMask(pthresh(i), NA, n, lambda, FRt); % generate Fourier mask
    Filt = Filt.*mask; % apply mask
    Irecon = ifftshift(ifft2(ifftshift(Filt))); % take iFFT
    LPFestimate = Irecon; % low pass filtered estimate
    Irecon = abs(Irecon); % positive/real values only 
    Irecon(isnan(Irecon)) = replvalue; % set nulls to zero
    Irecon(index) = FirstImage(index); % substitute back in sampled pixel values
    previous_image = Irecon; % most recent estimate
    
end
clkti = toc

% calculate SSIM score
cssim = ssim(normalize2D(Irecon), NOriginal); 

% look at data


figure(5); clf; tiledlayout(2,1,'tilespacing', 'compact');
nexttile; imagesc(VeryFirstImage); axis image; colormap(map); colorbar(); 
title('Original Image'); xticks([]); yticks([]);
nexttile; imagesc(normalize2D(Irecon)); axis image; colormap(map); colorbar(); 
title({'Reconstructed Image',sprintf('Fill-Factor: %0.1f%%, Cycles: %0.0f, Time: %0.3fs, SSIM: %0.4f', [ffarray{imind}, cyclenumber, clkti, cssim])});xticks([]); yticks([]);
% figure(6); clf; plot(sx, pthresh); axis tight; title('Threshold Function'); xlabel('Cycle'); ylabel('Threshold')
% figure(6); clf; plot(sx, maxthresh.*ones(length(sx))); axis tight; title('Threshold Function'); xlabel('Cycle'); ylabel('Threshold')


%% View original image
figure(1); clf; 
imagesc(VeryFirstImage); axis image; 
colormap(gray); xticks(''); yticks('')

%% Functions
function mask = generateMask(mult, NA, n, lambda, FRt)
    kc = mult*NA/(n*lambda); 
    CTF = rectpuls((FRt)./kc); 
    BeamSpatFreq = 1; 
    mask = BeamSpatFreq.*CTF; 
end