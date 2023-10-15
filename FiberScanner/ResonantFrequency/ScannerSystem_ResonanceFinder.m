% Scanner System Resonance Finder - Code Compiled for Matlab 2020a

% Author: Cameron Coleal
% Date: 11/24/2020

%   Description: This script can be used to run multiple scans in sequence
%   at set intervals. The data from each scan is saved and compiled at the
%   end of the script for comparative assesment. The overall goal is to
%   determine around which frequencies the piezo-fiber sees resonance in
%   both the x and y directions (determined independently)

%   Connections: May need to be adjusted based on your current setup
%   Requires the use of a photodetector (PDA36A2) and a PSD (Position
%   Sensitive Detector - PDP90A)

%   ai1: Photodetector intensity signal
%   PSD connections (11/24/20)
%       ai4: PSD connection #2 - Yaxis Pin
%       ai3: PSD connection #1 - Xaxis Pin
%       ai2: PSD connection #3 - Sum Pin


%% Resonance Range
xrange = (600:1:1000); % Set frequency scan range: X channel
yrange = (700:0.5:900); % Set frequency scan range: Y channel

%% Initialize DAQ
d = daq('ni');

% add output channels - Drive one output channel at a time
ch1 = addoutput(d,'Dev1','ao0','Voltage');% Analog Output for X Waveform
% ch2 = addoutput(d,'Dev1','ao1','Voltage');% Analog Output for Y Waveform

% add input channels
% ch3 = addinput(d,'Dev1','ai1','Voltage'); % Intensity Signal
ch4 = addinput(d,'Dev1','ai2','Voltage'); % PSD Sum Pin Signal
ch5 = addinput(d,'Dev1','ai3','Voltage'); % PSD X Pin Signal
ch6 = addinput(d,'Dev1','ai4','Voltage'); % PSD Y Pin Signal

% Set voltage range
ch5.Range = [-5,5]; 
ch6.Range = [-5,5];
%% Note: 
%  Pick either X or Y channel to analyze and uncomment the associated
%  channel, you will need to restart the script to then anaylze the next
%  channel (use "clc" and "clear" in the command line. Input channels are
%  provided but ideally you should be watching the response curves on the
%  NI DAQ Express software with this script. 
amplitude = 0.5; 
% amplitude = [0.1:0.1:1]; 

%% amplitude = 1; 

%% Loop through X frequencies and then y frequencies
% 
%  Navigate to folder
for ii = 1:length(xrange)
    % Round X (keep y the same, cycle over x values)

    %  Set DAQ parameters - xind
    fx = xrange(ii); 
    disp(fx)
    
    % DAQ Sampling Rate
    d.Rate = 83333; % Hz
    sampleDuration = 1; % 1/abs(fx-fy)
    nSamps = sampleDuration*d.Rate; 
    dt = 1/d.Rate;
    t = ([1:nSamps]-1)*dt; 

    % Initialization of waveforms
    OutX = amplitude*cos(2*pi*fx*t);
    
    % Run Scan
    [Signal] = readwrite(d,[OutX']); 

    % The variable "Signal" should be sized (:,2) where each index (1-2) corresponds to
    % the order we wrote the input data (X pin, Y pin)

    % Store the data  - this will be saved in the current folder unless the location
    % is changed

    % Combine time vector with signal matrix'
    SignalMat = [Signal{:,:},t'];
    
    
    % File saving details
    date = '2022May03'; 
    axis = 'X'; 
    scan_freq = strcat(string(fx), 'Hz.mat'); 
    folder_name = 'TestResonance_1sec_1HzStep\'; % You should be in this folder
    file_name = join([date, string(amplitude), 'ResScan', axis, scan_freq], '_'); % can change folder name

    % Save file
    save(strcat(folder_name, file_name), 'SignalMat')

end
% %% Reinitialize DAQ & rerun script with alt. analog channel uncommented before running next cell
% %% Run through y frequencies
% %  Navigate to folder
% for ii = 1:length(yrange)
%     % Round X (keep y the same, cycle over x values)
% 
%     %  Set DAQ parameters - xind
%     fy = yrange(ii); 
%     disp(fy)
%     
%     % DAQ Sampling Rate
%     d.Rate = 83333; % Hz
%     sampleDuration = 1; % 1/abs(fx-fy)
%     nSamps = sampleDuration*d.Rate; 
%     dt = 1/d.Rate;
%     t = ([1:nSamps]-1)*dt; 
% 
%     % Initialization of waveforms
%     OutY = amplitude*cos(2*pi*fy*t+pi/2);
% 
%     % Run the DAQ
%     [Signal] = readwrite(d,[OutY']); 
% 
%     % The variable "Signal" should be sized (:,2) where each index (1-2) corresponds to
%     % the order we wrote the input data (X pin, Y pin)
% 
%     % Store the data  - this will be saved in the current folder unless the location
%     % is changed
% 
%     % Combine time vector with signal matrix
%     SignalMat = [Signal{:,:}, t'];
% 
%     % File saving details
%     date = '2022Apr15'; 
%     axis = 'Y'; 
%     scan_freq = strcat(string(fy), 'Hz.mat'); 
%     folder_name = 'TestResonance_1sec_1HzStep\'; % You should be in this folder
%     file_name = join([date, 'ResScan', axis, scan_freq], '_'); % can change folder name
% 
%     % Save file
%     save(strcat(folder_name, file_name), 'SignalMat')
% 
% end

%% To plot recorded data run the "FrequencyResponsePlots.m" script
%% with amplitude loop
%% X output

for amp=1:length(amplitude)
    %  Navigate to folder
    for ii = 1:length(xrange)
        % Round X (keep y the same, cycle over x values)

        %  Set DAQ parameters - xind
        fx = xrange(ii); 
        disp(fx)

        % DAQ Sampling Rate
        d.Rate = 83333; % Hz
        sampleDuration = 1; % 1/abs(fx-fy)
        nSamps = sampleDuration*d.Rate; 
        dt = 1/d.Rate;
        t = ([1:nSamps]-1)*dt; 

        % Initialization of waveforms
        OutX = amplitude(amp)*cos(2*pi*fx*t);

        % Run Scan
        [Signal] = readwrite(d,[OutX']); 

        % The variable "Signal" should be sized (:,2) where each index (1-2) corresponds to
        % the order we wrote the input data (X pin, Y pin)

        % Store the data  - this will be saved in the current folder unless the location
        % is changed

        % Combine time vector with signal matrix
        SignalMat = [Signal{:,:},t'];


        % File saving details
        date = '2022Apr16'; 
        axis = 'X'; 
        scan_freq = strcat(string(fx), 'Hz.mat'); 
        folder_name = 'TestResonance_1sec_0.5HzStep\'; % You should be in this folder
        file_name = join([date, string(amplitude(amp)), 'ResScan', axis, scan_freq], '_'); % can change folder name

        % Save file
        save(strcat(folder_name, file_name), 'SignalMat')

    end
end

%% Y Output

for amp=1:length(amplitude)
        %  Navigate to folder
        for ii = 1:length(yrange)
                % Round X (keep y the same, cycle over x values)

        %  Set DAQ parameters - xind
        fy = yrange(ii); 
        disp(fy)

        % DAQ Sampling Rate
        d.Rate = 83333; % Hz
        sampleDuration = 1; % 1/abs(fx-fy)
        nSamps = sampleDuration*d.Rate; 
        dt = 1/d.Rate;
        t = ([1:nSamps]-1)*dt; 

        % Initialization of waveforms
        OutY = amplitude(amp)*cos(2*pi*fy*t+pi/2);

        % Run the DAQ
        [Signal] = readwrite(d,[OutY']); 

        % The variable "Signal" should be sized (:,2) where each index (1-2) corresponds to
        % the order we wrote the input data (X pin, Y pin)

        % Store the data  - this will be saved in the current folder unless the location
        % is changed

        % Combine time vector with signal matrix
        SignalMat = [Signal{:,:}, t'];

        % The variable "Signal" should be sized (:,2) where each index (1-2) corresponds to
        % the order we wrote the input data (X pin, Y pin)

        % Store the data  - this will be saved in the current folder unless the location
        % is changed

        % Combine time vector with signal matrix
%         SignalMat = [Signal{:,:},t'];+


        % File saving details
        date = '2022Apr16'; 
        axis = 'Y'; 
        scan_freq = strcat(string(fy), 'Hz.mat'); 
        folder_name = 'TestResonance_1sec_0.5HzStep\'; % You should be in this folder
        file_name = join([date, string(amplitude(amp)), 'ResScan', axis, scan_freq], '_'); % can change folder name

        % Save file
        save(strcat(folder_name, file_name), 'SignalMat')

    end
end
