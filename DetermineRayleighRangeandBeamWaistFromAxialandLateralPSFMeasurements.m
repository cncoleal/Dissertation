%% Code to plot resolution curves and find line of best fit
%  Cameron Coleal 
%  Spring, 2023 (JOSAA Revisions)

% units
mm = 1e-3; 
um = 1e-6; 

%% Load data arrays
% Lateral resolution data (23 rows by 4 columns)
% Position of Micrometer (mm) | Signal - 5s Average (V) | Signal - 5s Average (V) | Signal - 5s Average (V) | 
lateral = load('Lateral.mat').Lateral; 
Lx = lateral(:,1); 
Ly = lateral(:,2); 

% Axial resolution data (30 rowx by 2 columns)
% Position of Micrometer (mm) | Signal (V) | 
axial = load('Axial.mat').axial; 
Ax = axial(:,1); 
Ay = axial(:,2); 


%% ---------PROCESS LATERAL DATA------- %%
%% Subtract off DC -- prep for using Curve Fitter Application
LatDC = min(Ly); 
LX = (Lx - min(Lx))/mm; 
LY = Ly - LatDC; 

figure(1) 
clf; 
plot(LX, LY, '.-', 'LineWidth', 2, 'MarkerSize', 20)
xlabel('Position (um)')
ylabel('Signal (V)')
title('Lateral PSF Raw Data: DC Subtracted')
axis tight

%% Use matlab Curve Fitter Application to load data
%  Fit to custom equation: y = f(x) || y = I0^2*w0*pi*(sqrt(2)/2)*erfc(((x-x0)*sqrt(8))/w0);

%% Create figure to export using fit parameters
I0 =       0.09381;
w0 =         4.235;
x0 =         6.72;

latPSF =@(x) I0^2 * w0 * pi * (sqrt(2)/2) * erfc(((x-x0)*sqrt(8))/w0); 

x = linspace(min(LX), max(LX), 1000); 

fsize = 20; 
lwidth = 5; 
msize = 40; 


figure(2)
clf; 
tiledlayout(1,1,'TileSpacing','compact')
plot(x, latPSF(x), 'r', 'LineWidth', lwidth, 'DisplayName', 'Error Fxn Fit')
hold on 
plot(LX, LY,'k.-', 'MarkerSize', msize, 'DisplayName', 'Raw Data')
grid on 
xlabel('Position (um)','FontSize', fsize)
ylabel('Signal (V)','FontSize', fsize)
tt = sgtitle('Lateral PSF: Gaussian Complementary Error Fit');
set(tt, 'Horizontalalignment','center','fontsize', fsize)
legend('Location','northeast','FontSize', fsize-4)
xlim([min(LX), max(LX)])
set(gca, 'FontSize', 16)



%% ---------PROCESS AXIAL DATA------- %%
%% Clip axial data to ignore spherical aberration artifacts
[~, xind] = min(abs(Ax - 13.02)); 

AyClip = Ay(1:xind); 
AxClip = Ax(1:xind); 

figure(2) 
clf; 
subplot(131)
plot(Lx, Ly, 'LineWidth', 2)
xlabel('Position')
ylabel('Signal (V)')
title('Lateral PSF Raw Data')

subplot(132)
plot(Ax, Ay, 'LineWidth', 2)
xlabel('Position')
ylabel('Signal (V)')
title('Axial PSF Raw Data')

subplot(133)
plot(AxClip, AyClip, 'LineWidth', 2)
xlabel('Position')
ylabel('Signal (V)')
title('Axial PSF Clipped Data')

%% Subtract off DC -- prep for using Curve Fitter Application
AyDC = min(AyClip); 
AY = AyClip - AyDC; 
AyDC = Ay-min(Ay); 

figure(3) 
clf; 
plot(AxClip, AY, 'LineWidth', 2)
xlabel('Position')
ylabel('Signal (V)')
title('Axial PSF Clipped Data: DC Subtracted')


%% Create figure to export using fit parameters
DC =     -0.4739;
I0 =       1.487;
z0 =       0.09256;
zr =     0.08061;

fsize = 20; 
lwidth = 5; 
msize = 40; 

deltaDC = min(AyClip)-min(Ay); 

axPSF =@(z) I0^2 * ( 1 ./ (1 + ((2*(z-z0).^2) ./ (zr^2)) + (((z-z0).^4) / (zr^4)) )  ) + DC;
z = linspace(min(Ax), max(Ax), 1000); 
z = z-min(Ax); 

figure(4)
clf
tiledlayout(1,1,'TileSpacing','compact')
plot(z, axPSF(z)+deltaDC, 'r', 'LineWidth', lwidth, 'DisplayName', 'Fxn Fit')
hold on 
plot(Ax-min(Ax), (Ay-min(AyClip))+deltaDC,'k.-', 'MarkerSize', msize, 'DisplayName', 'Raw Data')
grid on 
xlabel('Position (um)','FontSize', fsize)
ylabel('Signal (V)','FontSize', fsize)

tt = sgtitle('Axial PSF: Lorenztian Fit');
set(tt, 'Horizontalalignment','center','fontsize', fsize)
legend('Location','northeast','FontSize', fsize-4)
axis tight
set(gca, 'FontSize', 16)
ylim([0, max((Ay-min(AyClip))+deltaDC)+.1])



%% ---------------SANITY CHECK RESULTS------- %%
%% Estimated beam waist from rayleigh range
nm = 1e-9; 
lambda = 785*nm; 
Zr = 80.6*um;   % rayleigh range
n = 1; 
est_w0 = (sqrt((lambda*Zr)/(pi*n)))/um % estimated beam waist (um)

%% Estimated rayleigh range from beam waist
nm = 1e-9; 
um = 1e-6; 
lambda = 785*nm; 
W0 = 4.23*um;   % beam waist
n = 1; 
est_zr = (((W0^2)*pi*n)/lambda)/um % estimated rayleigh range (um)
