%% Sweep Rotation of Input Drive Data
%  Cameron Coleal
%  3/25/21
%  The purpose of this script is to rotation through 360 degrees and
%  identify what angle minimizes coupling for driving X and Y frequencies
%  (i.e. identifying the eigen axes of the fiber).

%% Goal: Find Eigen Axes of Fiber
%% Initialize DAQ
d = daq('ni');

% add output channels
ch1 = addoutput(d,'Dev1','ao0','Voltage');% Analog Output for X Waveform
ch2 = addoutput(d,'Dev1','ao1','Voltage');% Analog Output for Y Waveform

% add input channels
% ch3 = addinput(d,'Dev1','ai1','Voltage'); % Intensity Signal
ch4 = addinput(d,'Dev1','ai2','Voltage'); % PSD Sum Pin Signal
ch5 = addinput(d,'Dev1','ai3','Voltage'); % PSD X Pin Signal
ch6 = addinput(d,'Dev1','ai4','Voltage'); % PSD Y Pin Signal

%% Change channel voltage range
ch5.Range = [-5,5]; 
ch6.Range = [-5,5]; 

%% Set amplitude
amp = 0.75; 
%%  Set DAQ parameters\

% 3 Frames/second
fx = 804; %Hz
fy = 807; %Hz
fo = 3; %GCF, Hz
nx = fx/fo; %119
ny = fy/fo; %118


k = 2; %symmetric about the origin for k=2,6

phix = 0; %phase on x is zero
phiy = (k*pi)/(4*nx); %phase on y

% Trace scan pattern
% Frequency for driving sinusoidal waveforms
% DAQ Sampling Rate
d.Rate = 83333; % Hz % note: can be higher than 62500 since only using 3 input channels

% Parameters
% Sample Duration (s)
 sampleDuration = 1; % 1/abs(fx-fy)

% Number of Samples 
nSamps = sampleDuration*d.Rate; 

% Sample Interval (s)
dt = 1/d.Rate;

% Create vector to hold sample data (s)
t = ([1:nSamps+1]-1)*dt;

% Output waveforms
DriveX = amp*cos(2*pi*nx*fo*t);
DriveY= amp*cos(2*pi*ny*fo*t+phiy);
%% Rotate DriveX
% Image Rotate
angle_range = (0:1:180); % angles to rotate through


%% Run Scan For X Axis
for ii = 1:length(angle_range)
    % Round X (keep y the same, cycle over x values)
   
    % set angle
    phi = angle_range(ii);
    % print(string(phi))
    % Initialization of waveforms
    zeroArray =  0.*[1:1:length(DriveX)]; 

    data = [DriveX; zeroArray];
    rotMat = [cosd(phi) -sind(phi); sind(phi) cosd(phi)];
    rotData = rotMat*data; 

    OutX = rotData(1,:); 
    OutY = rotData(2,:);
 
    % Run Scan
    [Signal] = readwrite(d,[OutX', OutY']);
    
    % Combine time vector with signal matrix
    SignalMat = [Signal{:,:},t'];
    
    
    % File saving details
    date = '2022May11'; % update with current date
    amplitude = amp; 
    axis = 'X'; 
    scan_angle = strcat(string(phi), 'Degrees.mat'); 
    folder_name = 'DriveFrequencyRotation\2022May11\3FramesPerSecond\90V\'; % This is the folder you will save to
    file_name = join([date, 'RotationScan', axis, amplitude, scan_angle], '_'); % can change folder name

    % Save file
    save(strcat(folder_name, file_name), 'SignalMat')
    

end

pause(2)

% Run Scan For Y Axis
for ii = 1:length(angle_range)
    % Round X (keep y the same, cycle over x values)
   
    % set angle
    phi = angle_range(ii);
    % Initialization of waveforms
%     DriveY = 0.5*cos(2*pi*fy*t); % X channel waveform
    zeroArray =  0.*[1:1:length(DriveY)]; 

    data = [DriveY; zeroArray];
    rotMat = [cosd(phi) -sind(phi); sind(phi) cosd(phi)];
    rotData = rotMat*data; 

    OutX = rotData(1,:); 
    OutY = rotData(2,:);
 
    % Run Scan
    [Signal] = readwrite(d,[OutX', OutY']);
    
    % Combine time vector with signal matrix
    SignalMat = [Signal{:,:},t'];
    
    
    % File saving details
    date = '2022May11';   % update with current date
    axis = 'Y'; 
    amplitude = amp; 
    scan_angle = strcat(string(phi), 'Degrees.mat'); 
    folder_name = 'DriveFrequencyRotation\2022May11\3FramesPerSecond\90V\'; % This is the folder you will save to
    file_name = join([date, 'RotationScan', axis,amplitude, scan_angle], '_'); % can change folder name

    % Save file
    save(strcat(folder_name, file_name), 'SignalMat')
    

end
 
%% Can rerun script with other frequency pairs to examine the impact of frequency selection on eigenaxis angle

% 1F/s
% fx = 709; %Hz
% fy = 710; %Hz
% fo = 1; %GCF, H

% % % 3F/s
% fx = 711; %Hz
% fy = 708; %Hz
% fo = 3; %GCF, Hz

% 
% % 4F/s
% fx = 708; %Hz
% fy = 712; %Hz
% fo = 4; %GCF, Hz
% 

% 6F/s
% fx = 714; %Hz
% fy = 708; %Hz
% fo = 6; %GCF, Hz

% 8F/s
% fx = 808; 
% fy = 800; 
% fo = 8; 

