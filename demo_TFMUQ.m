clearvars;

% Sampler parameters
alp0 = 1/(1e4)^2; % In units of Pa
bet0 = 1/(0.1)^2; % In units of disp field. inf to skip beta hyper-parameter
bet0 = inf;
params = struct();
params.maxIt = 50; 
params.burnIt = 1;
params.plotT = 0;

% saveIts = [];  % Empty array to skip saving to file
saveIts = struct();
saveIts.freq = 5; 
saveIts.savePath = './Results/test.mat';

% Gel parameters for TFM
E = 8e3; nu = 0.45; 
h = 150; % Height of the gel in microns


%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Read PIV with UQ %%%%%%%%%%%%%%%%%%%%%%%%%%
% Need to supply X,Y,U,V, Ustd and Vstd of the same size on [X,Y] grid 
% Ustd and Vstd are std. dev. estimates from UQ code
% U and V cannot contain nan or inf

calXY = 0.1628; % px -> um conversion
load('./Data/PIV.mat'); 
X = double(PIV.X * calXY); Y = double(PIV.Y * calXY); 
Ustd = double(PIV.Ustd); Vstd = double(PIV.Vstd);
U = double(PIV.U); V = double(PIV.V);

Ustd(isnan(Ustd))=max(Ustd(:)); Vstd(isnan(Vstd))=max(Vstd(:));
Ustd(Ustd==0)=max(Ustd(:)); Vstd(Vstd==0)=max(Vstd(:));
SigPIV = sparse([diag(Ustd(:)).^2 zeros(numel(Ustd)); zeros(numel(Ustd)) diag(Vstd(:)).^2]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UU = [U(:); V(:)]; % Reshape to vector

% Generate M matrix
if ~exist('M','var')
    tic;[M,~] = TfmMatrices2DFiniteThickLinButler(X,Y,X,Y,E,nu,1,1,1);toc;
end

% Smoothness prior
G = numgrid('S',size(X,1)+2);
D = delsq(G);
L = sparse([D zeros(size(D)); zeros(size(D)) D]); % Prior precision matrix

% % Global force balance prior (Uncomment the following line ) 
% L = L + ones(size(L)); 


%% Run the sampler 
[alpV,betV,TT,it] = MCMCGibbs_TFMUQ(M,UU,alp0,bet0,L,SigPIV,params,[],saveIts);


%% Visualization (TBD)


