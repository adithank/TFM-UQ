clearvars;
restoredefaultpath;

% Read the required images into a 3D array
imR = imread('./Data/Ref_10T12_pos9.tif');
imS = imread('./Data/Sess_10T12_pos9.tif');
im = cat(3,imR,imS);

PIV = PIVUQ('PIVConfig.json'); % Create a PIV object with the link to configuration 'json' file. Please refer to the json file for setting up PIV parameters
PIV.IMAGES = im; % Images to analyze
PIV.refFrames = 1; % Indices in the third dimension for reference frame
PIV.sessFrames = 2; % Indices in the third dimension for session frames (e.g. [2 3 4 ...])
PIV = PIV.runPIVMPUQ; % Run multi-pass PIV-UQ

% Post-process bootstrap samples using clustering
PIV = PIV.postProcess;

% Save results
out = PIV.output % ouput is a struct (see readme)
save ./Results/PIV.mat out;

%% Plotting example
plots.plotEgPIV(PIV.output(1)); 
