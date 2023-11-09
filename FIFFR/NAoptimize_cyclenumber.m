% Optimize a
function [cyclenumber_sampledDiff, optimized_cyclenumber_value] = optimize_cyclenumber(a, sigshift, maxthresh, minthresh, cyclenumber_list, Irecon, index, NOriginal, replvalue,  NA, n, lambda, FRt, simtype)
    % a_list: list of variable to test for sigmoid slope, a, value
    % sigshift: the value for shifting the sigmoid
    % maxthresh: the value for the maximum threshold
    % minthresh: the value for the minimum threshold
    % cyclenumber: the number of cycles to optimize for
    % Irecon: the normalized Irecon (with no NaN values!)
    % index: the logical array associated with sampled pixels
    % NOriginal: the original image, normalized
    % replvalue: replacement value for nulls
    % simtype: either 'psnr' or 'ssim' depending on which metric you want
    % to use

    % Variables to track: normalized
    normalize2D =@(M) (M-min(min(M)))/(max(max(M))-min(min(M))); 
    normalizeM = @(M) (M-min(M))/(max(M)-min(M)); 
    sigmoid =@(a, x, c, maxS, startS) (maxS-startS).*normalizeM(1./(1 + exp(-a.*(x-c))))+startS;  

    FirstImage = Irecon;

    SampledDiff = {}; 
    AbsSampledDiff = {}; 
    
    switch simtype
        case 'psnr'   
            for ii =1:length(cyclenumber_list)
                sx = linspace(-10, 10,cyclenumber_list(ii)); 
                pthresh = sigmoid(a, sx, sigshift, maxthresh, minthresh); 
            
                for i = 1:length(pthresh)
                    FTimage = fftshift(fft2(fftshift(Irecon))); 
                    Filt = FTimage; 
                    [mask] = generateMask(pthresh(i), NA, n, lambda, FRt); 
                    Filt = Filt.*mask;
                    Irecon = ifftshift(ifft2(ifftshift(Filt))); 

                    Irecon = abs(Irecon);
                    Irecon(isnan(Irecon)) = replvalue; % set nulls to zero
            
                    SampledDiff{i} =psnr(normalize2D(Irecon), NOriginal);
            
                    Irecon(index) = FirstImage(index); 

                    
                end
                carray = cell2mat(SampledDiff); 
                AbsSampledDiff{ii} = max(carray); 
            end


            cyclenumber_sampledDiff = cell2mat(AbsSampledDiff); 
            [~, cyclenumberind] = max(cyclenumber_sampledDiff); 
            optimized_cyclenumber_value = cyclenumber_list(cyclenumberind); 
    
        case 'ssim'
    
            for ii =1:length(cyclenumber_list)
                sx = linspace(-10, 10,cyclenumber_list(ii)); 
                pthresh = sigmoid(a, sx, sigshift, maxthresh, minthresh); 
            
                for i = 1:length(pthresh)
                    FTimage = fftshift(fft2(fftshift(Irecon))); 
                    Filt = FTimage; 
                    [mask] = generateMask(pthresh(i), NA, n, lambda, FRt); 
                    Filt = Filt.*mask;
                    Irecon = ifftshift(ifft2(ifftshift(Filt))); 

                    Irecon = abs(Irecon);
                    Irecon(isnan(Irecon)) = replvalue; % set nulls to zero
            
                    
                    Irecon(index) = FirstImage(index); 

                    SampledDiff{i} =ssim(normalize2D(Irecon), NOriginal);
                     
                end
                carray = cell2mat(SampledDiff); 
%                 AbsSampledDiff{ii} = min(1-carray); 
                AbsSampledDiff{ii} = max(carray);
            end

            cyclenumber_sampledDiff = cell2mat(AbsSampledDiff); 
            [~, cyclenumberind] = max(cyclenumber_sampledDiff); 
            optimized_cyclenumber_value = cyclenumber_list(cyclenumberind); 

    end


end