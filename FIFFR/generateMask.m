function mask = generateMask(mult, NA, n, lambda, FRt)
    kc = mult*NA/(n*lambda); 
    CTF = rectpuls((FRt)./kc); 
    BeamSpatFreq = 1; 
    mask = BeamSpatFreq.*CTF; 
end