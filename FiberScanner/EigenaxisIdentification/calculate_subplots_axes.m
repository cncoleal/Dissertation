function [logicalIm, corrxSig, corrySig, store, minX1, minY1,minX2, minY2,majX1, majY1,majX2, majY2] = calculate_subplots_axes(Signal, inputRot, viewPlots)
    sumpin = Signal{1,1}.SignalMat(:,1);
    xSig = 10*Signal{1,1}.SignalMat(:,2)./(2*sumpin);
    ySig = 10*Signal{1,1}.SignalMat(:,3)./(2*sumpin);
    
    % Create logical array
    V = ones(1, length(xSig)); 
    x = xSig; 
    y = ySig; 
    
    xg = linspace(min(x), max(x), 512); 
    yg = linspace(min(y), max(y), 512); 
    [Xg, Yg] = meshgrid(xg, yg); 
    Vg = griddata(x, y, V, Xg, Yg); 
    
%     figure(1)
%     imagesc(Vg)
    
    Vg(isnan(Vg))= 0; 
    logicalIm = Vg; 
    ctr = regionprops(logicalIm, 'centroid'); 
    theta = regionprops(logicalIm, 'orientation'); % ranges from -90 to 90
    major = regionprops(logicalIm, 'MajorAxisLength'); 
    minor = regionprops(logicalIm, 'MinorAxisLength'); 
    bb = regionprops(logicalIm, 'BoundingBox');
    
    store = cat(1,ctr.Centroid);
    majlen = cat(1,major.MajorAxisLength); 
    minlen = cat(1,minor.MinorAxisLength);
    angle = cat(1, theta.Orientation);
    bound = cat(1, bb.BoundingBox);

    % half the length of each axis
    halfMaj = majlen/2; 
    halfMin = minlen/2; 
    deg = angle(1); % dont take abs
    
    ctrx = store(1,1); 
    ctry = store(1,2); 
    
    % draw major axis
    delx = halfMaj*cosd(deg+90); 
    dely = halfMaj*sind(deg+90); 
    majX1 = ctrx + delx; 
    majX2 = ctrx - delx; 
    majY1 = ctry + dely; 
    majY2 = ctry - dely;
    
    % points to draw line of minor axis
    delyMin = halfMin*cosd(deg+90); 
    delxMin = halfMin*sind(deg+90); 
    minX1 = ctrx + (-1)*delxMin; 
    minX2 = ctrx - (-1)*delxMin; 
    minY1 = ctry + delyMin; 
    minY2 = ctry - delyMin; 
    
    % bounding box coordinates
    topleftx = bound(1,1); 
    toplefty = bound(1,2); 
    width = bound(1,3); 
    height = bound(1,4);
    
    topright = [topleftx+width, toplefty]; 
    bottomright = [topleftx+width, toplefty+height]; 
    bottomleft = [topleftx, toplefty+height]; 
    
    % convert xSig/ySig; 
    shiftxSig = xSig-min(xSig);
    normxSig = shiftxSig./(max(shiftxSig)-min(shiftxSig)) - 0.5; 
    corrxSig = width.*normxSig + width/2; 
    shiftySig = ySig-min(ySig); 
    normySig = shiftySig./(max(shiftySig)-min(shiftySig)) - 0.5; 
    corrySig = height.*normySig+ height/2; 



    if viewPlots == 1
        figure()
        imagesc(logicalIm)
        title(sprintf('%d Degrees',[inputRot]))
        hold on
        plot(corrxSig, corrySig)
        hold on
        plot(store(1,1), store(1,2), 'r*')
        hold on
        plot(minX1, minY1, 'r*')
        hold on 
        plot(minX2, minY2, 'r*')
        hold on 
        plot(majX1, majY1, 'g*')
        hold on 
        plot(majX2, majY2, 'g*')
        hold on
        plot([minX1, minX2], [minY1, minY2], 'r')
        hold on
        plot([majX1, majX2], [majY1, majY2], 'g')
        hold on
        plot(topleftx, toplefty, 'g+')
        hold on
        plot(topright(1), topright(2), 'g+')
        hold on
        plot(bottomright(1), bottomright(2), 'g+')
        hold on
        plot(bottomleft(1), bottomleft(2), 'g+')
        colormap('gray')
        axis off
    end
    
end