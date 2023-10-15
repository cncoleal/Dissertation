% Script to process resonance sweep at different amplitudes
% Cameron Coleal 
% March 08, 2022
% Note: naming convention is completely backward. All files labelled X are
% actually with the Y output intialized, and all files labelled Y are with
% the X output initiallized. Dont forget to change on graph names. 


%% create file name cells
xrange = (700:1:900); % Range: 300 - 400Hz
yrange = (700:1:900); % Range: 300 - 400Hz
amplist = 0.5%[0.1:0.1:1]; % range of amplitudes

% create file names
file_name_listX = {};
file_name_listY = {}; 

for amp = 1:length(amplist)
    for ii = 1:length(xrange) 
        % File saving details
        fx = xrange(ii); 
        fy = yrange(ii); 

        % X name list
        date = '2022May03'; 
        axis = 'X'; 
        scan_freq = strcat(string(fx), 'Hz.mat'); 
        folder_name = 'TestResonance_1sec_1HzStep\'; % The folder this folder is located in should be open on left side of matlab
        file_name = join([date, string(amplist(amp)),'ResScan', axis, scan_freq], '_'); % can change folder name
%         file_name = join([date, string(amplist(amp)),'ResScan',scan_freq], '_'); % can change folder name
        file_name_listX{amp, ii} = file_name; 


        % Y name list
        date = '2022Apr16'; 
        axis = 'Y'; 
        scan_freq = strcat(string(fy), 'Hz.mat'); 
        folder_name = 'TestResonance_1sec_0.5HzStep\'; % The folder this folder is located in should be open on left side of matlab
        file_name = join([date, string(amplist(amp)), 'ResScan', axis, scan_freq], '_'); % can change folder name
%         file_name = join([date, string(amplist(amp)), 'ResScan',  scan_freq], '_'); % can change folder name
        file_name_listY{amp, ii} = file_name; 

    end
end

%% Load all X files
%%  Signal Mat Details
%   SignalMat(:,1) = x Signal
%   SignalMat(:,2) = y Signal
%   SignalMat(:,3) = Time Vector

%% Load as cell array instead of structure
% Navigate to folder first

% Here -- SWAPPING BACK X AND Y CHANNELS
sXcell = {}; 

for amp=1:length(amplist)
    for ii=1:length(xrange)
       sXcell{amp, ii} = load(string(file_name_listX(amp, ii)));
    end
end

%% Read in Y data
sYcell = {}; 

for amp=1:length(amplist)
    for ii=1:length(xrange)
       sYcell{amp, ii} = load(string(file_name_listY(amp, ii)));
    end
end



%% Loop through and calculate Vpp
vppXx = {}; 
maxValXx = {}; 
vppXy = {}; 
maxValXy = {}; 
for amp=1:length(amplist)
    for ii=1:length(sXcell) 
        xdriveSUMread = sXcell{amp,ii}.SignalMat(:,1); 
        xdriveX = sXcell{amp,ii}.SignalMat(:,2);
        xdriveY = sXcell{amp,ii}.SignalMat(:,3);

        % apply sumpin
        xdriveXread = (10*xdriveX)./(2*xdriveSUMread); 
        xdriveYread = (10*xdriveY)./(2*xdriveSUMread); 

        maxXx = max(xdriveXread); 
        minXx = min(xdriveXread); 
        vppxx = abs(maxXx-minXx); 
        vppXx{amp, ii} = vppxx; 
        maxValXx{amp, ii} = maxXx; 

        maxXy = max(xdriveYread); 
        minXy = min(xdriveYread); 
        vppxy = abs(maxXy-minXy); 
        vppXy{amp, ii} = vppxy; 
        maxValXy{amp, ii} = maxXy;
    end
end

FreqRespXx = cell2mat(vppXx); 
MaxVoltXx = cell2mat(maxValXx); 
FreqRespXy = cell2mat(vppXy); 
MaxVoltXy = cell2mat(maxValXy); 

%% Plot Vpp X versus frequency
colors = {'r', 'b', 'g', 'k', 'm', '-.r', '-.b', '-.g', '-.k', '-.m'}; 

clf;
figure(1) 
for amp=1:length(amplist)
    voltval = amplist(amp)*120; 
    subplot(211)
    plot(xrange, FreqRespXx(amp,:), colors{1,amp},'DisplayName', strcat(string(voltval), 'V'))
    title('Drive X: X Frequency Response', 'FontSize', 24)
    xlabel('Driving Frequency', 'FontSize', 20)
    ylabel('Vpp over 1s (V)', 'FontSize', 20)
%     xlim([680, 740])
    set(gca, 'FontSize',18)
    grid minor
    legend()
    hold on
    
    subplot(212)
    plot(xrange, FreqRespXy(amp,:), colors{1,amp}, 'DisplayName', strcat(string(voltval), 'V'))
    title('Drive X: Y Frequency Response', 'FontSize', 24)
    xlabel('Driving Frequency', 'FontSize', 20)
    ylabel('Vpp over 1s (V)', 'FontSize', 20)
%     xlim([680, 740])
    grid minor
    legend()
    hold on
    set(gca, 'FontSize',18)
end

%% Plot Vpp X versus frequency --eigenaxis
colors = {'r', 'b', 'g', 'k', 'm', '-.r', '-.b', '-.g', '-.k', '-.m'}; 

clf;
figure(1) 
for amp=1:length(amplist)
    voltval = amplist(amp)*120; 
    subplot(211)
    plot(xrange, FreqRespXx(amp,:), colors{1,amp},'DisplayName', strcat(string(voltval), 'V'))
    title('Drive Piezotube X Axis: PSD X Axis Frequency Response', 'FontSize', 24)
    xlabel('Driving Frequency', 'FontSize', 20)
    ylabel('Vpp over 1s (V)', 'FontSize', 20)
%     xlim([680, 740])
    set(gca, 'FontSize',18)
    grid minor
    legend()
    hold on
    
    subplot(212)
    plot(xrange, FreqRespXy(amp,:), colors{1,amp}, 'DisplayName', strcat(string(voltval), 'V'))
    title('Drive Piezotube X Axis: PSD Y Axis Frequency Response', 'FontSize', 24)
    xlabel('Driving Frequency', 'FontSize', 20)
    ylabel('Vpp over 1s (V)', 'FontSize', 20)
%     xlim([680, 740])
    grid minor
    legend()
    hold on
    set(gca, 'FontSize',18)
end



%% Note for signal mat
% signalmat(:,1) - sumpin
% signalmat(:,2) - X
% signalmat(:,3) - Y
% signalmat(:,4) - time

%% Loop through and calculate Vpp
vppYy = {}; 
maxValYy = {}; 
vppYx = {}; 
maxValYx = {};  

for amp=1:length(amplist)
    for ii=1:length(sYcell) 
        ydriveSUMread = sYcell{amp,ii}.SignalMat(:,1); 
        ydriveX = sYcell{amp,ii}.SignalMat(:,2); 
        ydriveY = sYcell{amp,ii}.SignalMat(:,3); 
        
        % apply sumpin
        ydriveXread = (10*ydriveX)./(2*ydriveSUMread); 
        ydriveYread = (10*ydriveY)./(2*ydriveSUMread); 
        

        maxYy = max(ydriveYread); 
        minYy = min(ydriveYread); 
        vppyy = maxYy-minYy; 
        vppYy{amp, ii} = vppyy; 
        maxValYy{amp, ii} = maxYy; 

        maxYx = max(ydriveXread); 
        minYx = min(ydriveXread); 
        vppyx = maxYx-minYx; 
        vppYx{amp, ii} = vppyx; 
        maxValYx{amp, ii} = maxYx; 
    
    end
end

FreqRespYy = cell2mat(vppYy); 
MaxVoltYy = cell2mat(maxValYy); 
FreqRespYx = cell2mat(vppYx); 
MaxVoltYx = cell2mat(maxValYx); 

%% Plot Vpp X versus frequency
colors = {'r', 'b', 'g', 'k', 'm', '-.r', '-.b', '-.g', '-.k','-.m'}; 

clf; 
figure(2) 

for amp=1:length(amplist)
    voltval = amplist(amp)*120; 
    subplot(211)
    plot(yrange, FreqRespYx(amp,:), colors{1,amp},'DisplayName', strcat(string(voltval), 'V'))
    title('Drive Piezotube Y Axis: PSD X Axis Frequency Response', 'FontSize', 24)
    xlabel('Driving Frequency', 'FontSize', 20)
    ylabel('Vpp over 1s (V)', 'FontSize', 20)
%     xlim([680, 740])
    set(gca, 'FontSize',18)
    grid minor
    legend()
    hold on
    

    subplot(212)
    plot(yrange, FreqRespYy(amp,:), colors{1,amp},'DisplayName', strcat(string(voltval), 'V'))
    title('Drive Piezotube Y Axis: PSD Y Axis Frequency Response', 'FontSize', 24)
    xlabel('Driving Frequency', 'FontSize', 20)
    ylabel('Vpp over 1s (V)', 'FontSize', 20)
%     xlim([680,740])
    grid minor
    legend()
    hold on
    set(gca, 'FontSize',18)
    
end

%% Plot Vpp X versus frequency --eigenaxis
colors = {'r', 'b', 'g', 'k', 'm', '-.r', '-.b', '-.g', '-.k','-.m'}; 

clf; 
figure(2) 

for amp=1:length(amplist)
    voltval = amplist(amp)*120; 
    subplot(211)
    plot(yrange, FreqRespYx(amp,:), colors{1,amp},'DisplayName', strcat(string(voltval), 'V'))
    title('Drive Y: X Frequency Response', 'FontSize', 24)
    xlabel('Driving Frequency', 'FontSize', 20)
    ylabel('Vpp over 1s (V)', 'FontSize', 20)
%     xlim([680, 740])
    set(gca, 'FontSize',18)
    grid minor
    legend()
    hold on
    

    subplot(212)
    plot(yrange, FreqRespYy(amp,:), colors{1,amp},'DisplayName', strcat(string(voltval), 'V'))
    title('Drive Y: Y Frequency Response', 'FontSize', 24)
    xlabel('Driving Frequency', 'FontSize', 20)
    ylabel('Vpp over 1s (V)', 'FontSize', 20)
%     xlim([680,740])
    grid minor
    legend()
    hold on
    set(gca, 'FontSize',18)
    
end
%%
figure(4)
subplot(211)
plot(sXcell{4,150}.SignalMat(:,3), sXcell{4,150}.SignalMat(:,1))
title('Resonant Term')
subplot(212)
plot(sXcell{4,150}.SignalMat(:,3), sXcell{4,150}.SignalMat(:,2))
title('Cross Term')

figure(5)
subplot(211)
plot(sYcell{4,150}.SignalMat(:,3), sYcell{4,150}.SignalMat(:,1))
title('Cross Term')
subplot(212)
plot(sYcell{4,150}.SignalMat(:,3), sYcell{4,150}.SignalMat(:,2))
title('Resonant Term')

%% Plot frequency corresponding to the max value for each amplitude trace
MaxFreqXx = {}; 
XxInd={}; 
MaxFreqXy = {}; 
XyInd={}; 

MaxFreqYx = {}; 
YxInd={}; 
MaxFreqYy = {}; 
YyInd={}; 

for i=1:length(amplist)
    [x,y] = find(FreqRespXx(i,:)==max(FreqRespXx(i,:)));
    XxInd{i} = y; 
    MaxFreqXx{i} = xrange(y); 
    
    [x,y] = find(FreqRespXy(i,:)==max(FreqRespXy(i,:)));
    XyInd{i} = y;
    MaxFreqXy{i}  = xrange(y); 
    
    [x,y] = find(FreqRespYx(i,:)==max(FreqRespYx(i,:)));
    YxInd{i} = y;
    MaxFreqYx{i}  = yrange(y); 

    [x,y] = find(FreqRespYy(i,:)==max(FreqRespYy(i,:)));
    YyInd{i} = y;
    MaxFreqYy{i}  = yrange(y); 
end

%% view
linesize = 2; 
figure(10)

subplot(221)
plot(amplist*120, cell2mat(MaxFreqXx),'b*-', 'LineWidth', linesize)
title('MaxFreqXx vs Drive Voltage')
xlabel('Input Drive Voltage')
ylabel('Max Frequency')
xticks(amplist*120)
grid on 

subplot(223)
plot(amplist*120, cell2mat(MaxFreqXy),'b*-', 'LineWidth', linesize)
title('MaxFreqXy vs Drive Voltage')
xlabel('Input Drive Voltage')
ylabel('Max Frequency')
xticks(amplist*120)
grid on 

subplot(222)
plot(amplist*120, cell2mat(MaxFreqYx),'b*-', 'LineWidth', linesize)
title('MaxFreqYx vs Drive Voltage')
xlabel('Input Drive Voltage')
ylabel('Max Frequency')
xticks(amplist*120)
grid on 

subplot(224)
plot(amplist*120, cell2mat(MaxFreqYy),'b*-', 'LineWidth', linesize)
title('MaxFreqYy vs Drive Voltage')
xlabel('Input Drive Voltage')
ylabel('Max Frequency')
xticks(amplist*120)
grid on 

%% Show combined frequency response X/Y
Xcombined = zeros(size(sXcell)); 
Ycombined = zeros(size(sXcell)); 
for i=1:length(amplist)
    Xcombined(i,:) = sqrt(FreqRespXx(i,:).^2 + FreqRespXy(i,:).^2); 
    Ycombined(i,:) = sqrt(FreqRespYx(i,:).^2 + FreqRespYy(i,:).^2); 
end
%% Plot
c = ['r', 'b', 'g', 'k', 'm', '-.r', '-.b', '-.g', '-.k','-.m']; 
figure(11)
clf;
for amp=1:length(amplist)
    voltval = amplist(amp)*120;
    subplot(211)
    plot(xrange, Xcombined(amp, :), colors{1, amp}, 'DisplayName', strcat(string(voltval), 'V'))
    title('Drive Piezotube X Axis: Frequency Response', 'FontSize', 24)
    xlabel('Driving Frequency', 'FontSize', 20)
    ylabel('Vpp over 1s (V)', 'FontSize', 20)
    set(gca, 'FontSize',18)
    grid minor
    legend()
    hold on
    set(gca, 'FontSize',18)
    subplot(212)
    plot(yrange, Ycombined(amp, :), colors{1, amp},'DisplayName', strcat(string(voltval), 'V'))
    title('Drive Piezotube Y Axis: Frequency Response', 'FontSize', 24)
    xlabel('Driving Frequency', 'FontSize', 20)
    ylabel('Vpp over 1s (V)', 'FontSize', 20)
    set(gca, 'FontSize',18)
    grid minor
    legend()
    hold on
    set(gca, 'FontSize',18)
end

%% Calculate combined resonant frequencies
%% Plot frequency corresponding to the max value for each amplitude trace
MaxFreqX = {}; 
MaxFreqY = {}; 


for i=1:length(amplist)
    [x,y] = find(Xcombined(i,:)==max(Xcombined(i,:)));
    MaxFreqX{i} = xrange(y); 
    
    [x,y] = find(Ycombined(i,:)==max(Ycombined(i,:)));
    MaxFreqY{i}  = yrange(y); 
end

%% view
linesize = 2; 
figure(12)
clf;
subplot(121)
plot(amplist*120, cell2mat(MaxFreqX),'b*-', 'LineWidth', linesize)
title({'Drive X Axis','Max Frequency Response vs Drive Voltage'}, 'FontSize', 28)
xlabel('Input Drive Voltage (V)', 'FontSize', 24)
ylabel('Max Frequency (Hz)', 'FontSize', 24)
set(gca, 'FontSize',22)
xticks(amplist*120)
yticks(min(cell2mat(MaxFreqX)):0.5:max(cell2mat(MaxFreqX)))
xlim([12, 120])
grid on 

subplot(122)
plot(amplist*120, cell2mat(MaxFreqY),'b*-', 'LineWidth', linesize)
title({'Drive Y Axis','Max Frequency Response vs Drive Voltage'}, 'FontSize', 28)
xlabel('Input Drive Voltage (V)', 'FontSize', 24)
ylabel('Max Frequency (Hz)', 'FontSize', 24)
xticks(amplist*120)
yticks(min(cell2mat(MaxFreqY)):0.5:max(cell2mat(MaxFreqY)))
set(gca, 'FontSize',22)
xlim([12, 120])
grid on 



%% create table to view
startind = find(xrange==795); 
endind = find(xrange==815); 

freqvalues = xrange(startind:endind); 

freqXtable = zeros(length(amplist)+1, length(freqvalues)); 
freqYtable =zeros(length(amplist)+1, length(freqvalues)); 
freqXtable(1,:) = freqvalues; 
freqYtable(1,:) = freqvalues; 
for freq=1:length(amplist)
    freq
    freqXtable(freq+1,:) = Xcombined(freq, startind:endind); 
    freqYtable(freq+1,:) = Ycombined(freq, startind:endind); 
end
    
%% convert to table
Xtable = transpose(freqXtable); 
Ytable = transpose(freqYtable); 

%% create table to view for individual responses
startind = find(xrange==794); 
endind = find(xrange==818); 

freqvalues = xrange(startind:endind); 

freqXxtable = zeros(length(amplist)+1, length(freqvalues)); 
freqXytable = zeros(length(amplist)+1, length(freqvalues));
freqYxtable =zeros(length(amplist)+1, length(freqvalues)); 
freqYytable =zeros(length(amplist)+1, length(freqvalues)); 

freqXxtable(1,:) = freqvalues; 
freqXytable(1,:) = freqvalues;
freqYxtable(1,:) = freqvalues; 
freqYytable(1,:) = freqvalues; 

for freq=1:length(amplist)
    freqXxtable(freq+1,:) = FreqRespXx(freq, startind:endind); 
    freqXytable(freq+1,:) = FreqRespXy(freq, startind:endind); 
    
    freqYxtable(freq+1,:) = FreqRespYx(freq, startind:endind); 
    freqYytable(freq+1,:) = FreqRespYy(freq, startind:endind); 
end

%% transpose
Xxtable = transpose(freqXxtable); 
Xytable = transpose(freqXytable); 

Yxtable = transpose(freqYxtable); 
Yytable = transpose(freqYytable); 
    
%% plot X vs Y trace for 12V and 120V for each graph
figure(13)
subplot(121)
ind = 1;
signal = sXcell{ind, XxInd{ind}}; 
xview = (10*signal.SignalMat(:,2))./(2*signal.SignalMat(:,1)); 
yview = (10*signal.SignalMat(:,3))./(2*signal.SignalMat(:,1)); 
plot(xview, yview)
title({'Trace of Fiber Displacement at Resonance','Drive X Axis: 12V'}, 'FontSize', 28)
xlabel('PSD X Position (V)', 'FontSize', 24)
ylabel('PSD Y Position (V)', 'FontSize', 24)
set(gca, 'FontSize',22)
axis equal
xshift = -0.15; 
xlim([-1.5+xshift, 1.5+xshift])
ylim([-2, 2])
subplot(122)
ind = 10;
signal = sXcell{ind, XxInd{ind}}; 
xview = (10*signal.SignalMat(:,2))./(2*signal.SignalMat(:,1)); 
yview = (10*signal.SignalMat(:,3))./(2*signal.SignalMat(:,1)); 
plot(xview, yview)
title({'Trace of Fiber Displacement at Resonance','Drive X Axis: 120V'}, 'FontSize', 28)
xlabel('PSD X Position (V)', 'FontSize', 24)
ylabel('PSD Y Position (V)', 'FontSize', 24)
set(gca, 'FontSize',22)
axis equal
xlim([-1.5+xshift, 1.5+xshift])
ylim([-2, 2])
    
%% Y axis
figure(14)
subplot(121)
ind = 1;
signal = sYcell{ind, YyInd{ind}}; 
xview = (10*signal.SignalMat(:,2))./(2*signal.SignalMat(:,1)); 
yview = (10*signal.SignalMat(:,3))./(2*signal.SignalMat(:,1)); 
plot(xview, yview)
title({'Trace of Fiber Displacement at Resonance','Drive Y Axis: 12V'}, 'FontSize', 28)
xlabel('PSD X Position (V)', 'FontSize', 24)
ylabel('PSD Y Position (V)', 'FontSize', 24)
set(gca, 'FontSize',22)
axis equal
xshift = -0.15; 
xlim([-1.5+xshift, 1.5+xshift])
ylim([-2, 2])
subplot(122)
ind = 10;
signal = sYcell{ind, YyInd{ind}};  
xview = (10*signal.SignalMat(:,2))./(2*signal.SignalMat(:,1)); 
yview = (10*signal.SignalMat(:,3))./(2*signal.SignalMat(:,1)); 
plot(xview, yview)
title({'Trace of Fiber Displacement at Resonance','Drive Y Axis: 120V'}, 'FontSize', 28)
xlabel('PSD X Position (V)', 'FontSize', 24)
ylabel('PSD Y Position (V)', 'FontSize', 24)
set(gca, 'FontSize',22)
axis equal
xlim([-1.5+xshift, 1.5+xshift])
ylim([-2, 2])
    

%% Plot all traces at resonance
XxviewTable = {}; 
XyviewTable = {};
figure(15)
sgtitle('Trace of Fiber Displacement at Resonance')
for i=1:length(amplist)
    subplot(2,5,i)
    ind = i;
    signal = sXcell{ind, XxInd{ind}}; 
    xview = (10*signal.SignalMat(:,2))./(2*signal.SignalMat(:,1)); 
    yview = (10*signal.SignalMat(:,3))./(2*signal.SignalMat(:,1)); 
    XxviewTable{i} = xview; 
    XyviewTable{i} = yview; 
    plot(xview, yview)
    title(sprintf('Drive X Axis: %0.0f V', amplist(i)*120), 'FontSize', 18)
    xlabel('PSD X Position (V)', 'FontSize', 14)
    ylabel('PSD Y Position (V)', 'FontSize', 14)
    set(gca, 'FontSize',12)
    axis equal
    xshift = -0.15; 
    xlim([-1.5+xshift, 1.5+xshift])
    ylim([-2, 2])
end

%% 

YxviewTable = {}; 
YyviewTable = {};

figure(16)
sgtitle('Trace of Fiber Displacement at Resonance')
for i=1:length(amplist)
    subplot(2,5,i)
    ind = i;
    signal = sYcell{ind, YyInd{ind}}; 
    xview = (10*signal.SignalMat(:,2))./(2*signal.SignalMat(:,1)); 
    yview = (10*signal.SignalMat(:,3))./(2*signal.SignalMat(:,1));
    YxviewTable{i} = xview; 
    YyviewTable{i} = yview; 
    plot(xview, yview)
    title(sprintf('Drive Y Axis: %0.0f V', amplist(i)*120), 'FontSize', 18)
    xlabel('PSD X Position (V)', 'FontSize', 14)
    ylabel('PSD Y Position (V)', 'FontSize', 14)
    set(gca, 'FontSize',12)
    axis equal
    xshift = -0.15; 
    xlim([-1.5+xshift, 1.5+xshift])
    ylim([-2, 2])
end

%% test plot
figure(17)
plot(XyviewTable{1,1}, XxviewTable{1,1})
