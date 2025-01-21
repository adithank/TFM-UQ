function [xvec,yvec,tvec,U,V,varargout] = PIVRic2024(IMAGES,cfg_data,reference_frames,vecind,X0,Y0,T0,U0,V0,C0)
% Ricardo Serrano & Adithan Kandasamy 2024, all rights reserved
% Multipass 2D PIV for CPU. This function is to be used as one pass and
% should be invoked consecutive times in main code. If it is the first pass,
% X0,Y0,T0,U0,V0,C0 should be empty arrays.
%
% Inputs:
% - IMAGES: 3D-array where IMAGES(:,:,i) is an image from the experiment.
%           It should contain all images (reference and deformed).
% - cfg_data: Structure containing parameters of PIV configuration.
%       + cfg_data.Deformation.wdw_size: Window size
%       + cfg_data.Deformation.wdw_spacing: Window spacing
%       + cfg_data.Deformation.SubPixelInterpolationMode 'gaussian' or 'polynomial'
%           method for subpixel displacement
% - reference_frames: 1D array with index position of reference frames
%           within IMAGES stack. I.e. IMAGES(:,:,reference_frames) returns
%           frames to be used as reference.
% - vecind: 1D array with index position of frames to be processed within
%           IMAGES stack. I.e. IMAGES(:,:,vecind) returns frames to be analyzed
%           with PIV.
% - [X0,Y0,T0,U0,V0]: Spatio-temporal (X0,Y0,T0) coordinates of the deformation
%           field (U0,V0) from previous pass. This contains the information
%           for each window shift. For the first pass, pass all these
%           as empty variables '[]'.
% - C0: Cross-correlation matrices from previous run. They can be from
%           image quality assessment.
%
% Outputs:
% - [X,Y,T,U,V]: Spatio-temporal (X,Y,T) coordinates of the deformation
%           field (U,V). These are not saved to disk. User can decide
%           saving format in the main code.
%   Optional outputs:
%   A1, A2, A3, A4, A5  Parameters of the Parabolic subpixel interpolation
%                   f = A1*X^2 + A2*X + A3*Y^2 + A4*Y + A5
%   Cout: Cross-correlation for each image pair
%   IMS: Image windows
%   REFS: Reference windows

N_FRAMES = length(vecind);
N_REFS = length(reference_frames);
[IM_HEIGHT,IM_WIDTH, ~] = size(IMAGES);
WDW_SZ = cfg_data.Deformation.wdw_size;                  % window size;
WDW_SPC = cfg_data.Deformation.wdw_spacing;              % window distance
numelementsy=floor((IM_HEIGHT-WDW_SZ)/WDW_SPC)+1;
numelementsx=floor((IM_WIDTH-WDW_SZ)/WDW_SPC)+1;
N_WDWS = numelementsx*numelementsy;

% EDIT 1 : Precision of variables from json file
if isfield(cfg_data.Performance,'Precision')
    dataPrecision = cfg_data.Performance.Precision; % Used to be single
else
    dataPrecision = 'single';
end

if nargin < 4
    error('Not enough arguments');
elseif nargin == 4  % For first pass, 
    [X0,Y0,T0,U0,V0,C0] = deal([]);
elseif nargin <=9
    assert(nargin == 9);
    C0 = [];    
end

sx = 1:WDW_SPC:IM_WIDTH-WDW_SZ+1;
sy = 1:WDW_SPC:IM_HEIGHT-WDW_SZ+1;
st = vecind;
[xxx,yyy,ttt] = meshgrid(sx,sy,st);

xxx=xxx(:); yyy = yyy(:); ttt = ttt(:);

% Interpolate map to fit resolution of next pass
if ~isempty(X0)
    % Remove nans from deformation field of previous pass
    out = isnan(X0) | isnan(Y0) | isnan(T0) | isnan(U0) | isnan(V0);
    X0 = X0(~out); Y0 = Y0(~out); T0 = T0(~out); U0 = U0(~out); V0 = V0(~out);
    if N_FRAMES ==1
        Ump = griddata(double(X0 - WDW_SZ/2 + 1),double(Y0 - WDW_SZ/2 + 1),double(U0),xxx,yyy); % X0,Y0 are the centers of window whereas xini,yini are top left corners
        Vmp = griddata(double(X0 - WDW_SZ/2 + 1),double(Y0 - WDW_SZ/2 + 1),double(V0),xxx,yyy);
    else
        Ump = griddata(double(X0 - WDW_SZ/2 + 1),double(Y0 - WDW_SZ/2 + 1),double(T0*WDW_SPC),double(U0),xxx,yyy,ttt*WDW_SPC);
        Vmp = griddata(double(X0 - WDW_SZ/2 + 1),double(Y0 - WDW_SZ/2 + 1),double(T0*WDW_SPC),double(V0),xxx,yyy,ttt*WDW_SPC);
    end
else
    Ump = 0*xxx;
    Vmp = 0*xxx;
end


% At the edges, griddata produces nans (because at the edges it would need
% to extrapolate). Those nans are set to zero. Consider extrapolation of
% U0,V0 at edges.
Ump(isnan(Ump)|isnan(Vmp)) = 0;
Vmp(isnan(Ump)|isnan(Vmp)) = 0;

%
% Memory allocation helps performance
%
U = nan(N_WDWS*N_FRAMES,1,dataPrecision);
V = nan(N_WDWS*N_FRAMES,1,dataPrecision);
% Only needed for polynomial
if strcmp(cfg_data.Deformation.SubPixelInterpolationMode,'polynomial')
    A1 = nan(N_WDWS*N_FRAMES,1,dataPrecision);
    A2 = nan(N_WDWS*N_FRAMES,1,dataPrecision);
    A3 = nan(N_WDWS*N_FRAMES,1,dataPrecision);
    A4 = nan(N_WDWS*N_FRAMES,1,dataPrecision);
    A5 = nan(N_WDWS*N_FRAMES,1,dataPrecision);
end

if nargout == 13
    Cout = nan(WDW_SZ,WDW_SZ,N_WDWS*N_FRAMES,1,'double');
    IMS = zeros(WDW_SZ,WDW_SZ,N_WDWS*N_FRAMES,1,class(IMAGES));
    REFS = zeros(WDW_SZ,WDW_SZ,N_WDWS*N_FRAMES,1,class(IMAGES));
    
end

% Reshape the vector of top-left corner locations into 3D eases usage
% later on to obtain 3D array of windows.
xini_ref = permute(xxx,[3,2,1]); % With shape of [1 1 NN]
yini_ref = permute(yyy,[3,2,1]);

% Multipass consists in shifting top-left corner of the window
% (xini_ref,yini_ref) by Ump,Vmp
xini = round(xini_ref + permute(Ump,[3,2,1]));
yini = round(yini_ref + permute(Vmp,[3,2,1]));
tini = permute(ttt,[3,2,1]);


% Check for windows that will go out of bounds (oob_wdws) and drop them
% from analysis.
oob_wdws = (xini+WDW_SZ-1)>IM_WIDTH | (yini+WDW_SZ-1)>IM_HEIGHT | (xini <1) | (yini<1);
dropped_wdws = find(oob_wdws);
kept_wdws = find(~oob_wdws); % In linear indices 
N_KEPT_WDWS = length(kept_wdws);

%
% Automatic Batch Size Selection
%

a = whos('IMAGES');
b = whos('U');
MEM = a.bytes + 5 * b.bytes;
dum = ones(WDW_SZ,WDW_SZ,dataPrecision);
d = whos('dum');

% if ~isunix
%     t = memory;
%     MAX_RAM = min(t.MemAvailableAllArrays*0.9,cfg_data.Performance.MaxRAM*1e9);
% else
MAX_RAM = cfg_data.Performance.MaxRAM*1e9;
% end

if cfg_data.Performance.UseGPU
    MAX_GPU_MEM = cfg_data.Performance.MaxGPUMem*1e9;
end
BATCH_LENGTH_cpu = round(MAX_RAM/(20*d.bytes)); % Heuristic memory requirements
if cfg_data.Performance.UseGPU
    BATCH_LENGTH_gpu = round((MAX_GPU_MEM)/(20*d.bytes)); % Heuristic memory requirements
    BATCH_LENGTH = min([BATCH_LENGTH_cpu,BATCH_LENGTH_gpu]);
else
    BATCH_LENGTH = BATCH_LENGTH_cpu;
end
N_BATCHES = ceil(N_KEPT_WDWS/BATCH_LENGTH);

% Create the index matrices
[I,J] = meshgrid(double(0:WDW_SZ-1),double(0:WDW_SZ-1)); % (Y,X) coords for a single frame

%
%   Loop through batches of window pairs
% N_BATCHES
for ibatch = 1 : N_BATCHES
%     ibatch
    if ibatch ~= N_BATCHES
        curr_wdws_indx = [1:BATCH_LENGTH]'+(ibatch-1)*BATCH_LENGTH;
    else
        curr_wdws_indx = [1+(ibatch-1)*BATCH_LENGTH : N_KEPT_WDWS]';
        BATCH_LENGTH = length(curr_wdws_indx);
    end
    curr_wdws = kept_wdws(curr_wdws_indx); % Index of current windows in the batch
    
    Ises = repmat(I,[1 1 BATCH_LENGTH]); % [1:WDW_SZ,1:WDW_SZ,NBATCH] size
    Jses = repmat(J,[1 1 BATCH_LENGTH]);
        

    % 3D array of indexes to reshape IMAGES frame sequence into
    % interrogation windows of the batch.
    
    % yini(curr_wdws) is of shape [1 1 NBATCH] with NW corner of session
    % window coords. Matlab automatically performs repmat(a,[WDW_SZ,WDW_SZ,1]) while adding 
    
    INDimg = (Jses + yini(curr_wdws)) + (Ises + xini(curr_wdws)-1)*IM_HEIGHT + (tini(curr_wdws)-1)*IM_HEIGHT*IM_WIDTH;
    
    % Reshape images to array of windows, in GPU memory (IMAGES should be gpuArray before calling) or RAM.
    D = IMAGES(INDimg);
    
    %
    % Cross correlation (in Fourier space for perfomance)
    %
    D = fft2(D);
    D(1,1,:) = 0; % Set average = 0
    C = zeros(size(D),class(D));
    % Loop to accommodate Phase Average PIV (i.e. multiple reference frames)
    for iref = 1:N_REFS
        
        % 2D array of indexes to extract reference windows from IMAGES stack
        % Ises = Iref in this implementation
        INDref = (Jses + yini_ref(curr_wdws)) + (Ises + xini_ref(curr_wdws)-1)*IM_HEIGHT + (reference_frames(iref)-1)*IM_HEIGHT*IM_WIDTH;
        
        % NOTE: IMAGES should be gpuArray before calling for GPU computations
        R = IMAGES(INDref);
        
        R = fft2(R);
        R(1,1,:) = 0; % Set average = 0
        C = C + D.*conj(R); % Cross-correlation in Fourier space
    end
    % Bring back to real space
    C = fftshift(fftshift(real(ifft2(C)),1),2);
    
    % Handle peaks at zero
    %     Cint = 1/4*( C(1+WDW_SZ/2+1 ,1+WDW_SZ/2,:) + C(1+WDW_SZ/2,1+WDW_SZ/2+1,:) + ...
    %         C(1+WDW_SZ/2-1,1+WDW_SZ/2,:)+ C(1+WDW_SZ/2,1+WDW_SZ/2-1,:));
    %     C(1+WDW_SZ/2,1+WDW_SZ/2,:) = Cint;
    
    if ~isempty(C0)
        C = C-0.7*C0(:,:,curr_wdws);
    end
    if nargout == 14
        Cout(:,:,curr_wdws) = C;
        IMS(:,:,curr_wdws) = IMAGES(INDimg);
        REFS(:,:,curr_wdws) = IMAGES(INDref);   
    end
    clear INDimg
    clear D
    clear R
    
    %
    % Find location of maximum correlation
    %
    % Get peak location in C (the cross-correlation matrix) for subpixelar
    % interpolation. If correct, the following should be true: C(IND) = dum
    [dum, jg] = max(C,[],1); % Get row of max in columns
    [dum, ig] = max(dum,[],2); % Get location of max in X direction
    pixig = ig(:); % Pixelar velocity in X
    IG = sub2ind([WDW_SZ,BATCH_LENGTH],ig(:),[1:BATCH_LENGTH]');
    pixjg = jg(IG); % Pixelar velocity in Y
    IND = sub2ind([WDW_SZ,WDW_SZ,BATCH_LENGTH],pixjg(:),pixig(:),[1:BATCH_LENGTH]');
    
    %
    % Subpixelar interpolation
    %
    switch cfg_data.Deformation.SubPixelInterpolationMode
        case 'gaussian'
            [subpix_x,subpix_y] = gaussian_subpix(IND,WDW_SZ,C);
        case 'polynomial'
            [subpix_x,subpix_y,a1,a2,a3,a4,a5] = polynomial_subpix(IND,WDW_SZ,C);
            A1(curr_wdws) = single(a1);
            A2(curr_wdws) = single(a2);
            A3(curr_wdws) = single(a3);
            A4(curr_wdws) = single(a4);
            A5(curr_wdws) = single(a5);
    end
    peakxg = pixig(:) - (WDW_SZ/2 + 1) + subpix_x(:);
    peakyg = pixjg(:) - (WDW_SZ/2 + 1) + subpix_y(:);
    % Save results in X,Y,T,U,V variables and move on to next batch
    % They go in entries curr_wdws, leaving dropped_wdws as nan
    dum = (pixig < 2) | (pixjg < 2) | (pixig > (WDW_SZ-1)) | (pixjg > (WDW_SZ-1)); % | (abs(peakxg) > (WDW_SZ/2)) | (abs(peakyg) > (WDW_SZ/2));
    peakxg(dum) = nan;
    peakyg(dum) = nan;
    U(curr_wdws) = single(peakxg);
    V(curr_wdws) = single(peakyg);
    
end
% Reshape PIV outputs into 4D format. (y,x,frame,UQrepetition)
xvec = squeeze(single(sx));
yvec = squeeze(single(sy));
tvec = squeeze(single(vecind));

% X = reshape(single(xini_ref + WDW_SZ/2 - 1),numelementsy,numelementsx,N_FRAMES,cfg_data.Deformation.Nreps);
% Y = reshape(single(yini_ref + WDW_SZ/2 - 1),numelementsy,numelementsx,N_FRAMES,cfg_data.Deformation.Nreps);
% T = reshape(single(tini),numelementsy,numelementsx,N_FRAMES,cfg_data.Deformation.Nreps);

U = reshape(U + round(Ump),numelementsy,numelementsx,N_FRAMES);
V = reshape(V + round(Vmp),numelementsy,numelementsx,N_FRAMES);
if nargout > 6
    if strcmp(cfg_data.Deformation.SubPixelInterpolationMode,'polynomial')
        varargout{1} = reshape(A1,numelementsy,numelementsx,N_FRAMES);
        varargout{2} = reshape(A2,numelementsy,numelementsx,N_FRAMES);
        varargout{3} = reshape(A3,numelementsy,numelementsx,N_FRAMES);
        varargout{4} = reshape(A4,numelementsy,numelementsx,N_FRAMES);
        varargout{5} = reshape(A5,numelementsy,numelementsx,N_FRAMES);
    else
        varargout{1} = '';
        varargout{2} = '';
        varargout{3} = '';
        varargout{4} = '';
        varargout{5} = '';
    end
    varargout{6} = reshape(Cout,WDW_SZ,WDW_SZ,numelementsy,numelementsx,N_FRAMES);
    varargout{7} = reshape(IMS,WDW_SZ,WDW_SZ,numelementsy,numelementsx,N_FRAMES);
    varargout{8} = reshape(REFS,WDW_SZ,WDW_SZ,numelementsy,numelementsx,N_FRAMES);
    
end
% % Free-up GPU
% if cfg_data.Performance.UseGPU
%     gpuDevice(1);
% end
end
function [subpix_x,subpix_y] = gaussian_subpix(IND,WDW_SZ,C)
%
% Subpixelar interpolation (Gaussian)
%
% Protect against out of bounds
oob = (IND - WDW_SZ) < 1 | (IND+WDW_SZ) > numel(C);
duma = nan(length(IND),1);
dumb = duma;
dumc = duma;
dumd = duma;
dume = duma;
IND(oob) = [];

% Protect against negative values
duma(~oob) = C(IND); %(pixj,pixi)
dumb(~oob) = C(IND-WDW_SZ); %(pixj,pixi - 1)
dumc(~oob) = C(IND+WDW_SZ); %(pixj,pixi + 1)
dumd(~oob) = C(IND-1); %(pixj - 1,pixi)
dume(~oob) = C(IND+1); %(pixj + 1,pixi)
negvalue=(duma<0 | dumb<0 | dumc<0 | dumd<0 | dume<0);
duma(negvalue) = nan;
dumb(negvalue) = nan;
dumc(negvalue) = nan;
dumd(negvalue) = nan;
dume(negvalue) = nan;

% Subpixelar interpolation
f0 = reallog(duma);
f1 = reallog(dumb);%f1 = log(c(pixj,pixi-1));
f2 = reallog(dumc);%f2 = log(c(pixj,pixi+1));
subpix_x = (f1-f2)./(2*f1-4*f0+2*f2);
f1 = reallog(dumd);%f1 = log(c(pixj,pixi-1));
f2 = reallog(dume);%f2 = log(c(pixj,pixi+1));
subpix_y = (f1-f2)./(2*f1-4*f0+2*f2);

clear duma dumb dumc dumd dume
end
function [subpix_x,subpix_y,a1,a2,a3,a4,a5] = polynomial_subpix(IND,WDW_SZ,C)
%
% Subpixelar interpolation (Polynomial: A1*x^2 + A2*x + A3*y^2 + A4*y + A5)
%
dum11 = nan(length(IND),1);
dum12 = dum11;
dum13 = dum11;
dum21 = dum11;
dum22 = dum11;
dum23 = dum11;
dum31 = dum11;
dum32 = dum11;
dum33 = dum11;

% Protect against out of bounds
oob = (IND - WDW_SZ - 1) < 1 | (IND + WDW_SZ + 1) > numel(C);
IND(oob) = [];

dum11(~oob) = C(IND - WDW_SZ - 1); %(pixj - 1, pixi - 1)
dum12(~oob) = C(IND - 1);          %(pixj - 1, pixi    )
dum13(~oob) = C(IND + WDW_SZ - 1); %(pixj - 1, pixi + 1)
dum21(~oob) = C(IND - WDW_SZ);     %(pixj    , pixi - 1)
dum22(~oob) = C(IND);              %(pixj    , pixi    )
dum23(~oob) = C(IND + WDW_SZ);     %(pixj    , pixi + 1)
dum31(~oob) = C(IND - WDW_SZ + 1); %(pixj + 1, pixi - 1)
dum32(~oob) = C(IND + 1);          %(pixj + 1, pixi    )
dum33(~oob) = C(IND + WDW_SZ + 1); %(pixj + 1, pixi + 1)

a1 = 1/6 * ( dum11 - 2*dum12 + dum13 + dum21 - 2*dum22 + dum23 + dum31 - 2*dum32 + dum33 );
a2 = 1/6 * (-dum11 + dum13 - dum21 + dum23 - dum31 + dum33);
a3 = 1/6 * (dum11 + dum12 + dum13 - 2*dum21 - 2*dum22 - 2*dum23 + dum31 + dum32 + dum33);
a4 = 1/6 * (-dum11 - dum12 - dum13 + dum31 + dum32 + dum33);
a5 = 1/18* (-2*dum11 + 4*dum12 -2*dum13 + 4*dum21 + 10*dum22 + 4*dum23 - 2*dum31 + 4* dum32 - 2*dum33);

% Subpixelar interpolation
subpix_x = - a2./(2.*a1);
subpix_y = - a4./(2.*a3);
clear dum11 dum12 dum13 dum21 dum22 dum23 dum31 dum32 dum33

end