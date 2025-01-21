clearvars;
restoredefaultpath;

imR = imread('./Data/Ref_10T12_pos9.tif');
imS = imread('./Data/Sess_10T12_pos9.tif');

PIV = PIVUQ('PIVConfig.json');
PIV.IMAGES = cat(3,imR,imS);
PIV.refFrames = 1;
PIV.sessFrames = 2;
tic;PIV = PIV.runPIVMPUQ;toc

% Post-process bootstrap samplesusing clustering
PIV = PIV.postProcess;
