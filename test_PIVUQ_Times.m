% Test bootstrap times 
% clearvars;

im1 = imread('./Data/Ref_10T12_pos9.tif');
im2 = imread('./Data/Sess_10T12_pos9.tif');

cfg_data = loadJsonConfig('PIVConfig.json');

[X0,Y0,T0,U0,V0] = deal([]);

t=tic;
[xvec,yvec,tvec,~,U,V] = PIV_multipass_UQ_updatedUQ(cat(3,im1,im2),cfg_data,1,2,X0,Y0,T0,U0,V0,[]);
inbuilt = toc(t)

cfg_data.Deformation.Nreps = 1;

nB = 25;


t=tic;

imS = repmat(im2,[1 1 nB]);
N_RND = numel(imS);


% Indices to drop in the image
indS = randi(N_RND, [N_RND,1]);
indS = setdiff([1:N_RND],indS);
imS(indS) = 0;

[xvec,yvec,tvec,~,U,V] = PIV_multipass_UQ_updatedUQ(cat(3,im1,imS),cfg_data,1,1+[1:nB],X0,Y0,T0,U0,V0,[]);
bulk = toc(t)

%%
t=tic; 

parfor ii = 1 : nB
    im22 = im2;
    N_RND = numel(im2);

    rngID = rng('shuffle'); % For reproducibility

    % Indices to drop in the image
    indS = randi(N_RND, [N_RND,1]);
    indS = setdiff([1:N_RND],indS);
    im22(indS) = 0;
    
    [xvec,yvec,tvec,~,U,V] = PIV_multipass_UQ_updatedUQ(cat(3,im1,im2),cfg_data,1,2,X0,Y0,T0,U0,V0,[]);
end


Parforloop = toc(t)

%% GPU time

t=tic;
cfg_data.Performance.UseGPU = 1;
cfg_data.Performance.MaxGPUMem = 6;

im1G = gpuArray(im1);
im2G = gpuArray(im2);


cfg_data.Deformation.Nreps = nB;

[xvec,yvec,tvec,~,U,V] = PIV_multipass_UQ_updatedUQ(cat(3,im1,im2),cfg_data,1,2,X0,Y0,T0,U0,V0,[]);
GPUinbuilt = toc(t)


t=tic;

im1G = gpuArray(im1);
im2G = gpuArray(im2);

imS = repmat(im2G,[1 1 nB]);
N_RND = numel(imS);

% Indices to drop in the image
indS = randi(N_RND, [N_RND,1]);
indS = setdiff([1:N_RND],indS);
imS(indS) = 0;

cfg_data.Deformation.Nreps = 1;

[xvec,yvec,tvec,~,U,V] = PIV_multipass_UQ_updatedUQ(cat(3,im1,imS),cfg_data,1,1+[1:nB],X0,Y0,T0,U0,V0,[]);
GPUbulk = toc(t)

% Results 
% Inbuilt -43 s 
% GPU Inbuilt - 31 s 
% Bulk - 23 s 
% Forloop - 30 s
% Parfor (4 cores, 4 threads) - 11 s 
% GPUPB Bulk - 13.8 s 
