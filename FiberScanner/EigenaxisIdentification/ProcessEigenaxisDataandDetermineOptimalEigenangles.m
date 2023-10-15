%% Create figures for eigenaxis data 
%  Cameron Coleal
%  4/11/22

%%  Load data files
%% create file name cells
angle_range = (0:1:180); 

file_name_listX = {};
file_name_listY = {}; 

for ii = 1:length(angle_range) 
    % File saving details
    phi = angle_range(ii); 
    % File saving details
    date = '2022May11'; 
    axis = 'X_0.1'; 
    scan_angle = strcat(string(phi), 'Degrees.mat'); 
    folder_name = '12V/'; % You should be in this folder
    file_name = join([date, 'RotationScan', axis, scan_angle], '_'); % can change folder name
    file_name_listX{ii} = file_name; 
    
    
     % File saving details
    date = '2022May11'; 
    axis = 'Y_0.1'; 
    scan_angle = strcat(string(phi), 'Degrees.mat'); 
    folder_name = '12V/';% You should be in this folder
    file_name = join([date, 'RotationScan', axis, scan_angle], '_');
    file_name_listY{ii} = file_name; 

end
%% X data
%% Load as cell array instead of structure
% Navigate to folder first
sXcell = {}; 

for ii=1:length(angle_range)
   sXcell{ii} = load(string(file_name_listX(ii)));
end

%% Y data
%% Load as cell array instead of structure
sYcell = {}; 

for ii=1:length(angle_range)
   sYcell{ii} = load(string(file_name_listY(ii)));
end

%% process sXcell data
majorlengthx = {}; 
minorlengthx = {}; 

for ind = [1:length(sXcell)]
    viewPlots = 0; 
    [majlen, minlen] = calculate_major_minor_axes_final(sXcell(ind), ind, viewPlots);
    majorlengthx{ind} = majlen; 
    minorlengthx{ind} = minlen; 
end


%% process sYcell data
majorlengthy = {}; 
minorlengthy = {}; 

for ind = [1:length(sYcell)]
    viewPlots = 0; 
    [majlen, minlen] = calculate_major_minor_axes_final(sYcell(ind), ind, viewPlots);
    majorlengthy{ind} = majlen; 
    minorlengthy{ind} = minlen; 
end

%% FIGURE: show eigenaxis length vs input rotation for both x and y data
xmajor = cell2mat(majorlengthx);
xminor = cell2mat(minorlengthx);
fsize= 18; 
msize = 12; 
lsize = 2; 


figure(1)
clf
plot(xmajor, 'b', 'Linewidth', lsize)
xlabel('Degree of Input Rotation', 'FontSize', fsize)
ylabel('Major Axis Length', 'FontSize', fsize)
yyaxis left
hold on 
plot(xminor, 'r', 'Linewidth', lsize)
yyaxis right
ylabel('Minor Axis Length', 'FontSize', fsize)
% title('Eigenaxis Length vs Input X Rotation', 'FontSize', fsize)
title('Eigenaxis Length X: 804Hz with 12V Applied')
legend('Major Axis', 'Minor Axis', 'FontSize', fsize, 'Linewidth', lsize)
set(gca, 'FontSize', fsize)
axis tight
%% Only Minor Axis
figure(3)
clf
plot(xminor, 'r', 'Linewidth', lsize)
xlabel('Degree of Input Rotation', 'FontSize', fsize)
ylabel('Minor Axis Length', 'FontSize', fsize)
title('Eigenaxis Length vs Input X Rotation', 'FontSize', fsize)
legend('Minor Axis', 'FontSize', fsize, 'Linewidth', lsize)
set(gca, 'FontSize', fsize)
ylim([0 550])
xlim([0 180])

%% Y data
ymajor = cell2mat(majorlengthy);
yminor = cell2mat(minorlengthy);


figure(2)
clf;
plot(ymajor, 'b', 'Linewidth', lsize)
xlabel('Degree of Input Rotation', 'FontSize', fsize)
ylabel('Major Axis Length', 'FontSize', fsize)
yyaxis left
hold on 
plot(yminor, 'r', 'Linewidth', lsize)
yyaxis right
ylabel('Minor Axis Length', 'FontSize', fsize)
% title('Eigenaxis Length vs Input Y Rotation', 'FontSize', fsize)
title('Eigenaxis Length Y: 807Hz with 12V Applied')
legend('Major Axis', 'Minor Axis', 'FontSize', fsize, 'Linewidth', lsize)
set(gca, 'FontSize', fsize)
axis tight

%% Minor only
figure(4)
clf;
plot(yminor, 'r', 'Linewidth', lsize)
xlabel('Degree of Input Rotation', 'FontSize', fsize)
ylabel('Minor Axis Length', 'FontSize', fsize)
title('Eigenaxis Length vs Input Y Rotation', 'FontSize', fsize)
legend('Minor Axis', 'FontSize', fsize, 'Linewidth', lsize)
set(gca, 'FontSize', fsize)
ylim([0 550])
xlim([0 180])
%% show normalized graphs
normxmaj = xmajor./(max(xmajor)); 
normxmin = xminor./(max(xminor)); 
normymaj = ymajor./(max(ymajor)); 
normymin = yminor./(max(yminor)); 
fsize= 18; 
msize = 12; 
lsize = 2; 

figure(3)
clf;
plot(normxmaj, 'b', 'Linewidth', lsize)
xlabel('Degree of Input Rotation', 'FontSize', fsize)
ylabel('Major Axis Length (Normalized)', 'FontSize', fsize)
yyaxis left
hold on 
plot(normxmin, 'r', 'Linewidth', lsize)
yyaxis right
ylabel('Minor Axis Length (Normalized)', 'FontSize', fsize)
title('Eigenaxis Length vs Input X Rotation', 'FontSize', fsize)
legend('Major Axis', 'Minor Axis', 'FontSize', fsize)
set(gca, 'FontSize', fsize)


figure(4)
clf;
plot(normymaj, 'b', 'Linewidth', lsize)
xlabel('Degree of Input Rotation', 'FontSize', fsize)
ylabel('Major Axis Length (Normalized)', 'FontSize', fsize)
yyaxis left
hold on 
plot(normymin, 'r', 'Linewidth', lsize)
yyaxis right
ylabel('Minor Axis Length (Normalized)', 'FontSize', fsize)
title('Eigenaxis Length vs Input Y Rotation', 'FontSize', fsize)
legend('Major Axis', 'Minor Axis', 'FontSize', fsize)
set(gca, 'FontSize', fsize)

%% table
angle1x = find(normxmin(40:80) == min(normxmin(40:80)));
angle2x = find(normxmin(120:160) == min(normxmin(120:160))); 
angle1y = find(normymin(40:80) == min(normymin(40:80))); 
angle2y = find(normymin(120:160) == min(normymin(120:160))); 

table(angle1x+40-1, angle2x+120-1, angle1y+40-1, angle2y+120-1, 'VariableNames', {'Eigenaxis 1: X', 'Eigenaxis 2: X', 'Eigenaxis 1: Y', 'Eigenaxis 2: Y'})

%% Plot figures for 136 to 146
figure(5)
for ind = [55, 65, 70]
    viewPlots = 1; 
    [majlen, minlen] = calculate_major_minor_axes_final(sYcell(ind), ind, viewPlots);

end

%% black colormap
blackMap = [zeros(256, 1), zeros(256, 1), zeros(256, 1)]
colormap(blackMap);
%%
fsize= 18; 
msize = 12; 
count = 1; 
figure(6)
clf;
sgtitle('X Drive Eigenaxes', 'FontSize', fsize)
for ind = [30, 50, 57, 65, 70, 84]
%     ind = [30, 50,  57, 64, 84, ...
%         114, 134, 141, 148, 168]
    viewPlots = 0; 
    [logicalIm, corrxSig, corrySig, store, minX1, minY1,minX2, minY2,majX1, majY1,majX2, majY2] = calculate_subplots_axes(sYcell(ind), ind, viewPlots);
%     subplot(2,5,count)
    subplot(2,3,count)
    % set background
    background = logicalIm; 
    background(background~=0) = 0; 
    imagesc(background)
    colormap(blackMap)
    title(sprintf('%d Degrees',[ind]), 'FontSize', fsize)
    hold on
    h(1) = plot(corrxSig, corrySig)
    hold on
    plot(store(1,1), store(1,2), 'r*', 'MarkerSize', msize)
    hold on
    h(2) = plot(minX1, minY1, 'r*', 'MarkerSize', msize)
    hold on 
    plot(minX2, minY2, 'r*', 'MarkerSize', msize)
    hold on 
    h(3) = plot(majX1, majY1, 'g*', 'MarkerSize', msize)
    hold on 
    plot(majX2, majY2, 'g*', 'MarkerSize', msize)
    hold on
    plot([minX1, minX2], [minY1, minY2], 'r', 'LineWidth', 2, 'MarkerSize', msize)
    hold on
    plot([majX1, majX2], [majY1, majY2], 'g', 'LineWidth', 2, 'MarkerSize', msize)
    ax=gca; 
    ax.FontSize = fsize;
    axis off
    axis image
    count = count+1; 
    legend(h,'Scan Pattern', 'Minor Axis', 'Major Axis')
end


%% binary background Y 
fsize= 18; 
msize = 12; 
count = 1; 
figure(7)
clf;
sgtitle('Y Drive Eigenaxes', 'FontSize', fsize)
for ind = [30, 50, 57, 65, 70, 84]
%     ind = [30, 50,  57, 64, 84, ...
%         114, 134, 141, 148, 168]
    viewPlots = 0; 
    [logicalIm, corrxSig, corrySig, store, minX1, minY1,minX2, minY2,majX1, majY1,majX2, majY2] = calculate_subplots_axes(sYcell(ind), ind, viewPlots);
%     subplot(2,5,count)
    subplot(2,3,count)
    % set background
%     background = logicalIm; 
%     background(background~=0) = 0; 
    imagesc(logicalIm)
    colormap('gray')
    title(sprintf('%d Degrees',[ind]), 'FontSize', fsize)
    hold on
    h(1) = plot(corrxSig, corrySig)
    hold on
    plot(store(1,1), store(1,2), 'r*', 'MarkerSize', msize)
    hold on
    h(2) = plot(minX1, minY1, 'r*', 'MarkerSize', msize)
    hold on 
    plot(minX2, minY2, 'r*', 'MarkerSize', msize)
    hold on 
    h(3) = plot(majX1, majY1, 'g*', 'MarkerSize', msize)
    hold on 
    plot(majX2, majY2, 'g*', 'MarkerSize', msize)
    hold on
    plot([minX1, minX2], [minY1, minY2], 'r', 'LineWidth', 2, 'MarkerSize', msize)
    hold on
    plot([majX1, majX2], [majY1, majY2], 'g', 'LineWidth', 2, 'MarkerSize', msize)
    ax=gca; 
    ax.FontSize = fsize;
    axis off
    axis image
    count = count+1; 
    legend(h,'Scan Pattern', 'Minor Axis', 'Major Axis')
end



%% black background Y

fsize= 18; 
msize = 12; 
count = 1; 
figure(8)
clf;
sgtitle('Y Drive Eigenaxes', 'FontSize', fsize)
for ind = [30, 50,  57, 64, 84, ...
        114, 134, 141, 148, 168]
    viewPlots = 0; 
    [logicalIm, corrxSig, corrySig, store, minX1, minY1,minX2, minY2,majX1, majY1,majX2, majY2] = calculate_subplots_axes(sYcell(ind), ind, viewPlots);
    subplot(2,5,count)
    % set background
    background = logicalIm; 
    background(background~=0) = 0; 
    imagesc(background)
    colormap(blackMap)
    title(sprintf('%d Degrees',[ind]), 'FontSize', fsize)
    hold on
    h(1) = plot(corrxSig, corrySig)
    hold on
    plot(store(1,1), store(1,2), 'r*', 'MarkerSize', msize)
    hold on
    h(2) = plot(minX1, minY1, 'r*', 'MarkerSize', msize)
    hold on 
    plot(minX2, minY2, 'r*', 'MarkerSize', msize)
    hold on 
    h(3) = plot(majX1, majY1, 'g*', 'MarkerSize', msize)
    hold on 
    plot(majX2, majY2, 'g*', 'MarkerSize', msize)
    hold on
    plot([minX1, minX2], [minY1, minY2], 'r', 'LineWidth', 2, 'MarkerSize', msize)
    hold on
    plot([majX1, majX2], [majY1, majY2], 'g', 'LineWidth', 2, 'MarkerSize', msize)
    ax=gca; 
    ax.FontSize = fsize;
    axis off
    axis image
    count = count+1; 
    legend(h,'Scan Pattern', 'Minor Axis', 'Major Axis')
end


%%
%% Show zero and 180 degrees

fsize= 18; 
msize = 12; 
count = 1; 
figure(9)
clf;
sgtitle('90 and 180 Degree Input Rotation', 'FontSize', fsize)
for ind = [90, 180]
    viewPlots = 0; 
    [logicalIm, corrxSig, corrySig, store, minX1, minY1,minX2, minY2,majX1, majY1,majX2, majY2] = calculate_subplots_axes(sYcell(ind), ind, viewPlots);
    subplot(1,2,count)
    % set background
    background = logicalIm; 
    background(background~=0) = 0; 
    imagesc(background)
    colormap(blackMap)
    title(sprintf('%d Degrees',[ind]), 'FontSize', fsize)
    hold on
    h(1) = plot(corrxSig, corrySig)
    hold on
    plot(store(1,1), store(1,2), 'r*', 'MarkerSize', msize)
    hold on
    h(2) = plot(minX1, minY1, 'r*', 'MarkerSize', msize)
    hold on 
    plot(minX2, minY2, 'r*', 'MarkerSize', msize)
    hold on 
    h(3) = plot(majX1, majY1, 'g*', 'MarkerSize', msize)
    hold on 
    plot(majX2, majY2, 'g*', 'MarkerSize', msize)
    hold on
    plot([minX1, minX2], [minY1, minY2], 'r', 'LineWidth', 2, 'MarkerSize', msize)
    hold on
    plot([majX1, majX2], [majY1, majY2], 'g', 'LineWidth', 2, 'MarkerSize', msize)
    ax=gca; 
    ax.FontSize = fsize;
    axis off
    axis image
    count = count+1; 
    legend(h,'Scan Pattern', 'Minor Axis', 'Major Axis')
end

