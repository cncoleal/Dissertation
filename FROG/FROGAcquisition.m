%% This code acquires a FROG trace 
%  Cameron Coleal & Jesse Wilson
%  Last edited: 1/7/23

%  FORMAT: This code creates three cells to be used by the
%  FROGReconstruction_PCGPA.m code
%           - spectrum_pos (delay stage position)
%           - spectrum_wl  (wavelength vector from the spectrometer)
%           - spectrum_int (intensity vector measured from the
%           spectrometer)
%  PURPOSE: The purpose of this code is to acquired a FROG trace utilizing
%  a Newport Motor and a SeaBreeze spectrometer. Note the .DLL files for
%  both the motor/spectrometer are assumed to already be downloaded. 

%% Setup motor
% Initialize Linear Motor
if (~exist('motor','var'))
    motor = NewportMotor();
    motor.Options % view current motor specs
end

%% Setup the spectrometer
seaBreezeDir='C:\Program Files\Ocean Optics\SeaBreeze';
loadlibrary(fullfile(seaBreezeDir,'Library\SeaBreeze.dll'),...
    fullfile(seaBreezeDir,'API\SeaBreezeWrapper.h'),...
    'addheader',fullfile(seaBreezeDir,'API\DLLDecl.h'))

libfunctionsview SeaBreeze

% check api version
[ret1, ret2] = calllib('SeaBreeze','seabreeze_get_api_version_string',blanks(255),255)

% open connected spectrometer
ret = calllib('SeaBreeze','seabreeze_open_spectrometer',0,0)

% grab wavelength vector length
n = calllib('SeaBreeze','seabreeze_get_formatted_spectrum_length',0,0)

% get the wavelength vector
[err,len,wl] = calllib('SeaBreeze','seabreeze_get_wavelengths',0,0,zeros(n,1),n); 

%%  set the integration time of the spectrometer in microseconds
[err] = calllib('SeaBreeze','seabreeze_set_integration_time_microsec',0,0,1e5)

%% Parameters
c = 3e8;  %  speed of light
mm = 1e-3; % millimeters
nm = 1e-9; % nanometers
THz = 1e12; % terahertz
ps = 1e-12; % picoseconds

%% Setup parameters for motor

% Time overlap - grab current motor position if set at time overlap
toverlap = motor.Position; 

% parameters
c = 0.299792458; % speed of light in mm/ps

% Set scan range 
stepsize = 0.001; % step size of motor
minpos = toverlap-0.5; % minimum motor delay
maxpos = toverlap+0.5; % maximum motor delay

z = [minpos:stepsize:maxpos]; % array of motor positions

%% Clear variables - useful for multiple traces
clear spectrum_int spectrum_pos spectrum_wl

%% Run a continuous delay scan and plot spectrum at each delay 
% Setup GUI
hStopButton = uicontrol('Style','togglebutton',...
    'String','Stop','Callback',@(src,evt) stop(ni) );

% Set to rate
numscans =1; 
id = 0; 

% plot dtaa
figure(1)
clf;
title('Spectrum')
xlabel('Wavelength')
ylabel('Spectrometer Counts')
id = 0;

motor.Velocity = 0.5; % use faster return velocity to get to starting position
motor.Position = z(1);  % Motor goes to starting position

while id<numscans
    for ii=1:length(z)
        global data

        % Loop through the number of delay scans, n:
        motor.Position = z(ii);
        pause(0.5)

        % grab spectrometer data
        [err,len,s_wl] = calllib('SeaBreeze','seabreeze_get_formatted_spectrum',0,0,zeros(n,1),n);
        plot(wl(3:end),s_wl(3:end));
        ylim([0,4000]); % set the y axis limits for viewing
        xlim([490.5,540.5]); % set the x axis limits for viewing
        drawnow();
        
        % generate data cells
        spectrum_pos{ii} = z(ii);
        spectrum_wl{ii} = wl(3:end); 
        spectrum_int{ii} = s_wl(3:end); 
        
        pause(0.1)
        disp(z(ii)) % display current motor position
        
        
        id = id+1; 
    end
end


%% View the acquired data
nm = 1e-9; % nanometers
c = 3e8; % m/s
position = cell2mat(spectrum_pos); % motor position vector
tau = ((2.*position*mm)./c)/ps; % time delay vector, units ps
wavelength = cell2mat(spectrum_wl); % wavelength vector
wl = wavelength(:,1)';
frequency = (c./(wl.*nm))/THz; % frequency vector, units THz
spectrum = cell2mat(spectrum_int); 

wave = linspace(min(wl), max(wl), length(wl)); 
freq = linspace(min(frequency), max(frequency), length(wl)); 

IMwave = interp1(wl,spectrum,wave); % wave (rows) vs pos (cols)
IMfreq = interp1(frequency,spectrum,freq); % wave (rows) vs pos (cols)
IM1wave = IMwave; % pos (cols) vs wave (rows)

figure(1)
subplot(221)
imagesc(position,wave,IMwave)
ylabel('Wavelength (nm)')
xlabel('Position (mm)')

subplot(222)
imagesc(position,freq,IMfreq)
ylabel('Frequency (THz)')
xlabel('Position (mm)')

subplot(223)
imagesc(tau,wave,IMwave)
ylabel('Wavelength (nm)')
xlabel('Tau (ps)')

subplot(224)
imagesc(tau,freq,IMfreq) 
ylabel('Frequency (THz)')
xlabel('Tau (ps)')

%% Save data
filename = strcat('AcquiredFROGTrace');
save(filename, 'spectrum_pos', 'spectrum_wl', 'spectrum_int'); 




%% Useful functionalities %%
% -------------------------%
%% Reset motor position to overlap
motor.Position = toverlap; 


%% For independently using the spectrometer
%  set integration time in microseconds
[err] = calllib('SeaBreeze','seabreeze_set_integration_time_microsec',0,0,1e5)
%% View spectrum live
[err,len,wl] = calllib('SeaBreeze','seabreeze_get_wavelengths',0,0,zeros(n,1),n);
figure(2)
clf; 
while(1)
    [err,len,s_wl] = calllib('SeaBreeze','seabreeze_get_formatted_spectrum',0,0,zeros(n,1),n);
    plot(wl(3:end), s_wl(3:end), 'LineWidth', 3);
    xlim([490.5,540.5]);
    set(gca,'FontSize',40)
    drawnow();    
    pause(0.1)
    
end

