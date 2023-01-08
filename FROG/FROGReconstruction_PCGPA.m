%% This code uses the PCGPA FROG reconstruction method to estimate dispersion on ultrafast pulses
%  Cameron Coleal 
%  Last edited: 1/7/23

%  FORMAT: This code accepts three cells: 
%           - spectrum_pos (delay stage position)
%           - spectrum_wl  (wavelength vector from the spectrometer)
%           - spectrum_int (intensity vector measured from the
%           spectrometer)
%  To generate these vectors the code for FROG Acquisition can be used.
%  This file is titled: FROGAcquisition and utilizes a Newport Motor and
%  SeaBreeze spectrometer. 
%  PURPOSE: The purpose of this code is to process acquired FROG data and
%  estimate the first, second and third order dispersion values on an
%  ultrafast pulse. 

%% Units
THz = 1e12; % terahertz
ps = 1e-12; % picoseconds
fs = 1e-15; % femtoseconds
nm = 1e-9;  % nanometers
um = 1e-6;  % micrometers
mm = 1e-3;  % millimeters

%% Load data and convert into matrices
spectrum_pos = cell2mat(spectrum_pos); % motor position
spectrum_wl = cell2mat(spectrum_wl);  % wavelength 
spectrum_int = cell2mat(spectrum_int); % spectrums

%% Create the FROG trace
wave_c = 400*nm; % Set the center wavelength of the SFG spectrum
tau_bw =3.8*ps; % Specify the time delay bandwidth 
N = 512; % Specify the number of pixels the FROG 2D array should have, NxN

% create the 2D FROG Array
[tau, freq, time, RealFrog] = CreateFROGConstraint(wave_c, tau_bw, N, spectrum_pos, spectrum_wl, spectrum_int); 

% Display the generated FROG trace
figure(1)
clf;
imagesc(tau./ps, freq./THz, RealFrog)
xlabel('Delay (ps)')
ylabel('Frequency (THz)')
title('Experimental Data')
colorbar()

%% Subtract off the background 

% Apply a median filter if the FROG trace is noisy
FrogMed = medfilt2(RealFrog); 

% Specify a pixel intensity corresponding to the highest background value
px_background = 185; 

% subtract the background
FF = FrogMed-px_background; 

% set all background pixels to zero
FF(FF<0) = 0; 

% normalize the FROG 
IFrog = (FF-min(min(FF)))./(max(max(FF))-min(min(FF))); 

% View the Original vs Background subtracted traces with a custom colormap
figure(2)
clf;
subplot(121)
imagesc(tau./ps, freq./THz,IFrog)
cmap = colormap();
cmap(1,:)=0;
colormap(cmap);
xlabel('Delay (ps)')
ylabel('Frequency (THz)')
title('Experimental Data with Background Removed')
colorbar()
subplot(122)
imagesc(tau./ps, freq./THz,IFrog)
cmap = colormap();
cmap(1,:)=0;
colormap(cmap);
xlabel('Delay (ps)')
ylabel('Frequency (THz)')
title('Experimental Data with Background Removed')
colorbar()
ylim([-20,20])


% View the Original vs Background subtracted traces with the default colormap
figure(3)
clf; 
subplot(121)
imagesc(tau./ps, freq./THz, sqrt(abs(IFrog)))
xlabel('Delay (ps)')
ylabel('Frequency (THz)')
title('Sqrt(Abs(Experimental Data with Background Removed))')
colormap default
colorbar()
subplot(122)
imagesc(tau./ps, freq./THz, sqrt(abs(IFrog)))
xlabel('Delay (ps)')
ylabel('Frequency (THz)')
title('Sqrt(Abs(Experimental Data with Background Removed))')
colormap default
colorbar()
ylim([-20,20])


%% Generate initial guess
% Time vector
% Gaussian 
E0 = 1; 
taup = 141*fs; % Estimated FWHM of Pulse
taug = taup/1.177; 
Et =@(t) E0*exp(-2*((t/taug).^2)); 

% start with simple gaussian estimates for probe
P = Et(time); 

figure(4)
clf;
plot(time/ps, P)
xlabel('Delay (ps)')
ylabel('Frequency (THz)')
title('Initial Probe Guess')

%% Run FROG Loop 
fL = 15; % the upper and lower frequency limits [THz] to display in the viewing windows while the reconstruction is running
tL = 1; % the upper and lower time delay limits [ps] to display in the viewing windows while the reconstruction is running
Probe = P/(max(P)); 

% Set an error threshold if you wish to run the reconstruction until an
% error is met
% error_threshold = 0.1; % UNCOMMENT
 
% Set an initial error value
% error = 100; % UNCOMMENT

% Keep track of iterations by using the iteration variable if wanted
% iteration = 1; % UNCOMMENT

% while error > error_threshold
for iteration=1:200
    % Box 2: first outer product
    Outer1 = Probe.'*Probe; 

    figure(5)
    clf;
    hStopButton = uicontrol('Style','togglebutton','String','Stop',...
        'Position',[500,50,40,20],'Callback',@stopCallback);
    subplot(4,4,[1 5])
    imagesc(tau/ps,time/ps,abs(Outer1))
    title('Box 2: OuterProduct')
    xlabel('Delay')
    ylabel('Time')
    axis square

    % Box 3: Row Rotate to the Time Domain
    % Go through each row and shift by 1
    shift = 0; 
    for row = 1:length(Probe) 
        Outer1(row, :) = circshift(Outer1(row, :), shift); 
        shift = shift-1; 
    end

    subplot(4,4,[2 6])
    imagesc(tau/ps,time/ps,abs(Outer1))
    title('Box 2: OuterProduct Shift')
    xlabel('Delay (ps)')
    ylabel('Time (ps)')
    axis square

     % Box 4: Take the fourier transform along columns
    FFT = fftshift(fft(fftshift(Outer1, 2),[],1),1);

    % view
    subplot(4,4,[3 7])
    imagesc(tau/ps,freq/THz,abs(FFT))
    title('Box 4: Estimated FROG')
    xlabel('Delay (ps)')
    ylabel('Frequency (THz)')
    axis square
    
    %%Box 5: Apply Intensity Constraint
    % Replace modulus of guess with actual modulus
    % phase array from guess
    phase = angle(FFT); 

    % calculate error
    error = sqrt(sum(sum((FFT-IFrog).^2)))/(N); 
 
    % New guess
    Guess = sqrt(abs(IFrog)).*exp(1i*phase); 
    
    subplot(4,4,[4 8])
    imagesc(tau/ps,freq/THz,abs(Guess))
    title('Box 5: Replaced Modulus - New Guess')
    xlabel('Delay (ps)')
    ylabel('Frequency (THz)')
    axis square

     % Box 6: Take the inverse fourier transform along rows
     iFFT = fftshift(ifft(ifftshift(Guess,1),[],1),2); 

    subplot(4,4,[9 13])
    imagesc(tau/ps,time/ps,abs(iFFT))
    title('Box 6: iFFT')
    xlabel('Delay (ps)')
    ylabel('Time (ps)')
    axis square
        
    % reverse circshift
    shift = 0; 
    for row = 1:length(Probe) 
        iFFT(row, :) = circshift(iFFT(row, :), shift); 
        shift = shift+1; 
    end
        
    subplot(4,4,[10 14])
    imagesc(tau/ps,time/ps,abs(iFFT))
    title('Box 7: iFFT Shift')
    xlabel('Delay (ps)')
    ylabel('Time (ps)')
    axis square
    
    % Box 8: SVD
    [U, W, V] = svd(iFFT, 0); 
    
    % new gate 
    Probe = transpose(U(:,1)); 

    Probe = Probe./max(Probe); 

    subplot(4,4,[11 15])
    plot(time/ps,abs(Probe),'k')
    yyaxis right
    plot(time/ps,unwrap(angle(Probe)))
    title({sprintf('Probe: Iteration %0.0f', iteration),sprintf('Error: %0.5f', error)})
    xlabel('Time (ps)')
    axis square
    
    subplot(4,4,12)
    imagesc(tau/ps,freq/THz,abs(FFT))
    title('Box 4: Zoomed In')
    xlabel('Delay (ps)')
    ylabel('Frequency (THz)')
    axis square
    ylim([-fL, fL])
    xlim([-tL, tL])
    
    subplot(4,4,16)
    imagesc(tau/ps,freq/THz,abs(Guess))
    title('Box 5: Zoomed In')
    xlabel('Delay (ps)')
    ylabel('Frequency (THz)')
    axis square
    ylim([-fL, fL])
    xlim([-tL, tL])
    
%     iteration = iteration + 1; % UNCOMMENT
 
end



%% Plot final frog estimate versus reconstruction
fL =20; % frequency axis limits for viewing
tL = 2; % time delay axis limits for viewing

figure(6)
subplot(121)
imagesc(tau/ps,freq/THz,abs(FFT).^2)
title('Final FROG Estimate')
xlabel('Delay (ps)')
ylabel('Frequency (THz)')
axis square
ylim([-fL, fL])
xlim([-tL, tL])

subplot(122)
imagesc(tau/ps,freq/THz,abs(IFrog))
title('Experimental Data')
xlabel('Delay (ps)')
ylabel('Frequency (THz)')
axis square
ylim([-fL, fL])
xlim([-tL, tL])
colormap default

%% Save Variables (if desired)
% filename = strcat('FROGLoopVariables'); % UNCOMMENT
% save(filename, 'error', 'FFT', 'freq', 'IFrog', 'iteration', 'Probe', 'RealFrog', 'tau', 'tau_bw', 'time', 'wave_c');  % UNCOMMENT

%% Take the Fourier Transform of the estimated Probe pulse
c = 3e8; % speed of light
freq_c = c/(wave_c); % center frequency  
ProbeF = fftshift(fft(fftshift(Probe))); % Fourier transform of the probe pulse

% create a new freq vector
Nv = length(time); 
dt = time(2)-time(1); 
dv = 1/Nv/dt; 
Fw = dv*Nv; 

% this is same frequency basis as freq vector
ff = linspace(-Fw/2, Fw/2, Nv);

% wavelength basis 
wvl = (c./(ff+freq_c))/nm; 
wv = linspace(min(wvl), max(wvl), length(ProbeF)); 

% View temporal and spectral estimates
figure(7)
clf; 
tiledlayout(1,2)
nexttile
plot(time/ps,abs(Probe),'k')
ylabel('Temporal Intensity (arb.)')
yyaxis right
plot(time/ps,unwrap(angle(Probe)),'r')
title('Temporal Pulse Profile')
xlabel('Time (ps)')
ylabel('Spectral Phase (rad)')
xlim([min(time/ps), max(time/ps)])

nexttile
plot(freq/THz,abs(ProbeF), 'k')
ylabel('Spectral Intensity (arb.)')
yyaxis right
plot(freq/THz, unwrap(angle(ProbeF)),'r')
xlabel('Frequency (THz)')
ylabel('Spectral Phase (rad)')
title('Spectral Pulse Profile')
xlim([min(freq/THz), max(freq/THz)])

%% Clip the data frame for fitting the phase function -- want to avoid edges of the spectrum
MatchPhase = unwrap(angle(ProbeF)); % phase
Amp = abs(ProbeF); % amplitude

%  Clip the spectrum
Fwid =2*THz; % frequency spectrum width
[~,indL] = min(abs(freq+Fwid)); 
[~,indH] = min(abs(freq-Fwid)); 

% non-normalized data
PhaseYd = MatchPhase(indL:indH); 
AmpYd = Amp(indL:indH); 

% frequency axis clipped
freqXd = freq(indL:indH)/THz; 

figure(8)
clf; 
plot(freqXd, AmpYd, 'k')
ylabel('Magnitude')
yyaxis right
plot(freqXd, PhaseYd, 'b.')
ylabel('Spectral Phase (rad)')
xlabel('Frequency (THz)')
xlim([min(freqXd), max(freqXd)])
title('Freq Domain: Probe')


%% Data Fit
% fit the phase with a 4th order polynomial 
order = 4; 

% Raw data fit
fitobject = fit(freqXd', PhaseYd', 'poly4');  

if order ==1
    % 1ST ORDER POLY FIT
    fitfun =@(x,FO) FO.p1*x+ FO.p2;
elseif order == 2
    % 2ND ORDER POLY FIT
    fitfun =@(x,FO) FO.p1*x.^2 + FO.p2*x+ FO.p3;
elseif order ==3
     % 2ND ORDER POLY FIT
     fitfun =@(x,FO) FO.p1*x.^3 + FO.p2*x.^2 + FO.p3*x + FO.p4;
elseif order == 4
    % 4TH ORDER POLY FIT
    fitfun =@(x,FO) FO.p1*x.^4 + FO.p2*x.^3 + FO.p3*x.^2 + FO.p4*x + FO.p5; 
end

% fit function
f = fitfun(freqXd, fitobject); 

% view data 
Pcenter = (max(PhaseYd)-min(PhaseYd))/2+min(PhaseYd);
f = (f-min(f))/(max(f)-min(f))*(max(PhaseYd)-min(PhaseYd)); 
fcenter = (max(f)-min(f))/2+min(f);
fshift =Pcenter-fcenter; 
f = f+fshift; 

figure(9)
clf; 
plot(freqXd, AmpYd, 'k','DisplayName', 'Probe Amplitude')
ylabel('Magnitude')
yyaxis right
plot(freqXd, PhaseYd, 'b.','DisplayName', 'Probe Phase')
hold on 
plot(freqXd, f, 'r','DisplayName', 'Probe Phase Fit')
ylabel('Spectral Phase (rad)')
xlabel('Frequency (THz)')
xlim([min(freqXd), max(freqXd)])
title({'Freq Domain: Probe',['Fxn: ',func2str(fitfun)]})
legend('Location','NorthWest')


% create table with values
% extract vars
if order ==1
    % 1ST ORDER POLY FIT
    pv = [fitobject.p2, fitobject.p1]; 
    rN = {'p2'; 'p1'}; 
elseif order == 2
    % 2ND ORDER POLY FIT
    pv = [fitobject.p3,fitobject.p2, fitobject.p1]; 
    rN = {'p3';'p2'; 'p1'};
elseif order ==3
     % 3RD ORDER POLY FIT
    pv = [fitobject.p4,fitobject.p3,fitobject.p2, fitobject.p1]; 
    rN = {'p4';'p3';'p2'; 'p1'};
elseif order == 4
    % 4TH ORDER POLY FIT
    pv = [fitobject.p5,fitobject.p4,fitobject.p3,fitobject.p2, fitobject.p1]; 
    rN = {'p5';'p4';'p3';'p2'; 'p1'};
end

% flip vectors
Pb = fliplr(pv)';

% Calculate GD
% ps2 to fs2 - unit transformations
psfs = 1000; % picosecond to femtosecond
ps2fs2 = psfs^2; % picoseconds^2 to femtoseconds^2
ps3fs3 = psfs^3; % picoseconds^3 to femtoseconds^3

% Group Delay
GD = fitobject.p4/(2*pi); % [ps]

% Group Delay Dispersion
GDD = fitobject.p3/(4*pi^2)*ps2fs2; % [fs^2]

% Third Order Dispersion
TOD = fitobject.p2/(8*pi^3)*ps3fs3; % [fs^3]

% Reformat table vectors
Pb = [Pb; GD; GDD; TOD]; 
el = length(rN); 
rowNam = flipud(rN); 
rowNam{el+1} = 'GD (ps)'; 
rowNam{el+2} ='GDD (fs^2)'; 
rowNam{el+3} = 'TOD (fs^3)'; 

% make table
T  = table(Pb, 'RowNames', rowNam, 'VariableNames', {'Values'}); 

% view table
figure(10)
clf; 
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);







%% This function takes the measured FROG data and generates a 2D FROG array usable by the reconstruction algorithm. 
function [tau, freq,time, IFrog] = CreateFROGConstraint(wave_c, tau_bw, N, spectrum_pos, spectrum_wl, spectrum_int)
    % INPUTS: 
    %       - wave_c: [scalar] center wavelength of the measured spectrum
    %       - tau_bw: [scalar] specified bandwidth to keep in the 2D FROG array
    %       - N: [integer] number of pixels in the length and width of the 2D FROG
    %       arrray (NxN)
    %       - spectrum_pos: [double] 1D array of the motor positions 
    %       - spectrum_wl: [double] 2D array of the wavelengths measured by
    %       the spectrometer at each delay
    %       - spectrum_in: [double] 2D array of the measured spectrums at
    %       each delay


    % OUTPUTS: 
    %       - tau: [double] 1D array of the temporal positions of the motor
    %       delay within the specified bandwidth, size 1xN, zero-centered. 
    %       - freq: [double] 1D array of the wavelengths measured by the
    %       spectrometer within the specified bandwidth, size 1xN, zero-centered
    %       - time: [double] 1D array, corresponding time axis to the freq
    %       vector, size 1xN, zero-centered. 


    % units
    mm = 1e-3; % millimeters
    nm = 1e-9; % nanometers
    ps = 1e-12; % picoseconds
    THz = 1e12; % terahertz

    % speed of light 
    c = 3e8; % [m/s]

    % delay position vector from motor
    position_raw = spectrum_pos; 
    position = linspace(min(position_raw), max(position_raw), length(position_raw)); % evenly spaced position delay vector

    % delay time vector from motor
    tau_raw = ((2.*position*mm)./c); % [s]
    tau = linspace(min(tau_raw), max(tau_raw), length(tau_raw)); % evenly spaced time delay vector

    % wavelength array from spectrometer
    wavelength = spectrum_wl; 

    % wavelength vector, reformated so 1D and evenly spaced
    wl_raw = wavelength(:,1)';
    wl = linspace(min(wl_raw), max(wl_raw), length(wl_raw)); % evenly spaced wavelength vector

    % frequency vector Hz, created from wavelength vector
    frequency_raw =(c./(wl.*nm)); % [Hz]
    frequency = linspace(min(frequency_raw), max(frequency_raw), length(frequency_raw)); % evenly spaced frequency vector

    % spectral intensity array, 
    spectrum = spectrum_int; 

    % Sampling spectrum array onto different bases

    % First sample raw data onto smooth wavelength basis
    IM0 = interp1(wl_raw,spectrum, wl,'linear',  0);
    % Second sample raw data onto a smooth position basis
    IM1 = (interp1(position_raw, IM0', position, 'linear', 0))'; 

    % Third interpolate onto a frequency axis from wavelength
    IM2 = interp1(frequency_raw,IM1, frequency,'pchip', 0);
    % Fourth interpolate onto a time delay axis from position delay
    IMAGE = (interp1(tau_raw,IM2', tau,'pchip',  0))';

    % Clip the specified bandwidths, then interpolate the data
    % center frequency
    freq_c = (c/(wave_c)); 

    % find indices corresponding to the center frequency
    [~,fmid] = min(abs(frequency-freq_c));

    % center delay 
    lineoutFreqMid = IMAGE(fmid, :); % grab the lineout at the center frequency
    xind = linspace(1, length(tau),length(tau)); 

    % center of mass of center lineout
    cm = sum(lineoutFreqMid.*xind)/sum(lineoutFreqMid); % find the center of mass
    tau_c = tau(1,round(cm));  % set the center delay position to the center of mass of the frequency lineout

    % tau clip - find the lower and upper times of the clipped time delay vector
    tau_ll = tau_c-(tau_bw/2); % lower limit
    tau_ul = tau_c+(tau_bw/2); % upper limit

    % find indices corresponding to the clipped bandwidth
    [~,tlower] = min(abs(tau-tau_ll)); % lower index
    [~,tmid] = min(abs(tau-tau_c)); % middle index
    [~,tupper]= min(abs(tau-tau_ul)); % upper index

    % check that specified bandwidth is supported by data collected
    if tau_bw > (max(tau)-min(tau))
        sprintf('Reduce Tau Bandwidth by %0.07f ps', ((max(tau)-min(tau))-tau_bw)/ps)
        sprintf('Tau Bandwidth Available is %0.02f ps', (max(tau)-min(tau))/ps)
    else
        disp('Tau Bandwidth Okay')
    end

    % CLIP the time delay axis
    TAUclip = tau(tlower:tupper); 
    TAU = linspace(min(TAUclip), max(TAUclip), N);

    dTau = TAU(2)-TAU(1); 
    dv = 1/N/dTau; 
    Fw = dv*N; 

    % find the upper and lower frequencies of the clipped frequency axis
    fLL = freq_c-(Fw/2); % lower frequency 
    fUL = freq_c+(Fw/2); % upper frequency

    % find limits
    [~,flower] = min(abs(frequency-fLL)); % lower index
    [~,fupper] = min(abs(frequency-fUL)); % upper index

    % clip the 2D Frog
    IMclip = IMAGE(flower:fupper, tlower:tupper); 

    % create the new frequency axis vector post-clipping
    FREQclip = frequency(flower:fupper); 

    % Interpolate the 2D frog onto the clipped tau/freq axes
    FREQ = linspace(fLL, fUL, N); 
    finalIM = interp1(FREQclip, IMclip, FREQ,'linear',  0); 
    FROGreal = (interp1(TAUclip, finalIM', TAU,'linear',  0))'; 

    % shift to center at zero 
    tau = TAU-tau_c; % centered delay
    freq = FREQ-freq_c; % centered frequency

    % generate time vector
    df = freq(2)-freq(1); 
    dt = 1/N/df; 
    Tw = dt*N; 
    time = linspace(-Tw/2, Tw/2, N); % time vector

    % Final FROG 2D array for reconstruction
    IFrog = FROGreal; 
    % X axis is delay
    % Y axis is frequency 
end

function stopCallback(src,event)
    stop()
end