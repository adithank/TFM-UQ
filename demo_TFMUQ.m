clearvars;

% MCMC sampler parameters
alp0 = 1/(1e4)^2; % In units of Pa (prior scale parameter, alpha)
bet0 = 1/(0.1)^2; % In units of disp field (epsilon_beta). inf to skip beta hyper-parameter
% bet0 = inf; 
params = struct();
params.maxIt = 25; % maximum MCMC iterations after burn-in
params.burnIt = 20; % Burn-in iterations to discard
params.plotT = 1; % Plot evert step 

% saveIts = [];  % Empty array to skip saving to file
saveIts = struct();
saveIts.freq = 5; 
saveIts.savePath = './Results/MCMC_save.mat';

% Gel parameters for TFM
E = 8e3; nu = 0.45; 
h = 150; % Height of the gel in microns


%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Read PIV with UQ %%%%%%%%%%%%%%%%%%%%%%%%%%
% Need to supply X,Y,U,V, Ustd and Vstd of the same size on [X,Y] grid 
% Ustd and Vstd are std. dev. estimates from UQ code
% U and V cannot contain nan or inf

calXY = 0.1628; % px -> um conversion
PIV=load('./Results/PIV.mat'); 
PIV=PIV.out;

[X,Y] = meshgrid(PIV.xvec, PIV.yvec);
X = double(X * calXY); Y = double(Y * calXY); 
Ustd = double(PIV.UStd)*calXY; Vstd = double(PIV.VStd) * calXY;
U = double(PIV.UPost) * calXY; V = double(PIV.VPost) * calXY;

Ustd(isnan(Ustd))=max(Ustd(:)); Vstd(isnan(Vstd))=max(Vstd(:)); % Remove "bad windows" from bootstrap analysis
Ustd(Ustd==0)=max(Ustd(:)); Vstd(Vstd==0)=max(Vstd(:)); % Remove 0's 
SigPIV = sparse([diag(Ustd(:)).^2 zeros(numel(Ustd)); zeros(numel(Ustd)) diag(Vstd(:)).^2]); % Create sigma_PIV covariance matrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UU = [U(:); V(:)]; % Reshape to vector

% Generate M matrix
if ~exist('M','var')
    tic;[M,~] = MoP.TfmMatrices2DFiniteThickLinButler(X,Y,X,Y,E,nu,1,1,1);toc;
end

% Smoothness prior
G = numgrid('S',size(X,1)+2);
D = delsq(G);
L = sparse([D zeros(size(D)); zeros(size(D)) D]); % Prior precision matrix

% % Global force balance prior (Uncomment the following line ) 
% L = L + ones(size(L)); 


%% Run the MCMC sampler 
[alpV,betV,TT,it] = MCMC.MCMCGibbs_TFMUQ(M,UU,alp0,bet0,L,SigPIV,params,[],saveIts);

%% Save results
TFMUQ = struct();
TFMUQ.X = X;
TFMUQ.Y = Y;
TFMUQ.alpV = alpV;
TFMUQ.betV = betV;
TTXUQ = reshape(TT(1:end/2,:),[size(X) size(TT,2)]);
TTYUQ = reshape(TT(1+end/2:end,:),[size(X) size(TT,2)]);
TFMUQ.TTXUQ = TTXUQ;
TFMUQ.TTYUQ = TTYUQ;
TFMUQ.it=it;
TFMUQ.params = params;
save('./Results/TFMUQ.mat','TFMUQ','U','V','Ustd','Vstd');
        
%% Visualization 
clf; plots.plotEgTFMUQ(TFMUQ); 

