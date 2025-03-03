%% Calculate 2D TFM matrix for finite thickness using DFT
% Lin et al, Mechanosesning of substrate thickness (Butler lab)
% https://journals.aps.org/pre/pdf/10.1103/PhysRevE.82.041918
%
% INPUT :
%   [X,Y] - Output grid points (u)
%   [Xc,Yc] - Measurement or traction grid points (t)
%   E - Young's modulus in the units of output
%   nu - Poisson's ratio
%   h - Thickness of the gel
%   realDom - 1 for Real (multiplied by DFT), 0 for Fourier space (FAST)
%   forwardOp - Need M matrix (t to u map) too? 
% OUTPUT : 
%   M - Full convolution matrix of the form u = Mt; 
%       u being displacement
%   Minv - M inverse so that t = M^-1 u
%   and t, traction forces. u :  m x 1 ; t : n x 1; M : m x n
% 
%


function [M,Minv] = TfmMatrices2DFiniteThickLinButler(XD,YD,XT,YT,E,nu,h,realDom,forwardOp,Zero0thmode)

    if nargin < 10
        Zero0thmode = 1;
    end
    dxD = unique(XD);
    NxD = length(dxD);
    dxD = dxD(2) - dxD(1);
    dyD = unique(YD);
    NyD = length(dyD);
    dyD = dyD(2) - dyD(1);
    dxT = unique(XT);
    NxT = length(dxT);
    dxT = dxT(2) - dxT(1);
    dyT = unique(YT);
    NyT = length(dyT);
    dyT = dyT(2) - dyT(1);
    
    LxD = XD(end) - XD(1); % Domain length of uniform grid
    LyD = YD(end) - YD(1);
    LxT = XT(end) - XT(1); % Domain length of uniform grid
    LyT = YT(end) - YT(1);
    
    % Checking for grid irregularities
    if ~( numel(XD) == NxD*NyD ) | ~(numel(YD) == NxD*NyD )
        error(' Check for redundancy or uniformness of displacement grid');
    end
    if ~( numel(XT) == NxT*NyT)  | ~( numel(YT) == NxT*NyT )
        error(' Check for redundancy or uniformness of traction grid');
    end
    % Since the Fourier TFM equation is for each mode, the Nx and Ny of
    % both disp and traction grid has to be the same
    if NxD ~= NxT | NyD ~= NyT
        error(' Number of grid points are not the same in both grids');
    end
    
    k = [-NxD/2:NxD/2-1]*((2*pi)/(dxD*NxD));
    l = [-NyD/2:NyD/2-1]*((2*pi)/(dyD*NyD));
    [K,L] = meshgrid(k,l);
    
    Ks = fftshift(K); % alpha
    Ls = fftshift(L); % beta
    K2s = (Ks.^2 + Ls.^2);
    q = sqrt(K2s); % q is the scalar, modulus of 2D wave vector (k,l)
    
    % Defintions
    c = cosh(q.*h);
    s = sinh(q.*h);
    
    % Constructing M matrix for vector u_hat = [u_k ; u_l];
    f1 = E.* (c.*q)./(2.*(1+nu).*s);
    f2 = ( E./(2.*(1-nu.^2).*q.*s) ) .* ...
            ( (3-4*nu)*nu.*s.*(c.^2) - (1-nu).*c.*q.*h + (1-2*nu).^2.*s + s.*(q.*h).^2 ) ./ ...
                    ( (3-4*nu).*s.*c + q.*h );
    %     MMinv = kron(f1 ,eye(NxD*2,NyD*2)) + kron( f2, [Ks.^2 Ks.*Ls; Ks.*Ls Ls.^2]);
    %     MMinv = kron( f1 ,eye(NxD*2,NyD*2)) + kron( f2.*Ks.^2 , repmat([1 0;0 0],NxD, NyD)) + ...
    %         kron( f2.*Ls.^2 , repmat([0 0;0 1],NxD, NyD)) + kron( f2.*Ks.*Ls , repmat([0 1;1 0],NxD, NyD) );
    Af = diag(reshape(f1 + f2.*Ks.^2,[],1));
    BCf = diag(reshape(f2.*Ks.*Ls,[],1));
    Df = diag(reshape(f1+f2.*Ls.^2,[],1));
    MMinv = [Af BCf; BCf Df];
    clear Af BCf Df;
    MMinv(isnan(MMinv)) = 0;
    if Zero0thmode
        MMinv(1,:) = 0;
        MMinv(1+end/2,:) = 0;
    end
        
    % TRY ONLY
%     f1(1) = f1(2);
%     f2(1) = f2(2);
    
    g1 = 1./f1;
    g2 = -f2./(f1.*(f1+q.^2 .*f2));
    
    Ag = diag(reshape(g1 + g2.*Ks.^2,[],1));
    BCg = diag(reshape(g2.*Ks.*Ls,[],1));
    Dg = diag(reshape(g1+g2.*Ls.^2,[],1));
    MM = [Ag BCg; BCg Dg];
    MM(isnan(MM)) = 0;
    if Zero0thmode
        MM(1,:)=0;
        MM(1+end/2,:) = 0;
    end
    clear Ag BCg Dg;
    
    % Real space
    
    if realDom == 0
        if forwardOp == 1
            M = MM;
        end
        Minv = MMinv;
    elseif realDom == 1
        [Fv,Finv] = MoP.genDFTmatricesFor2DVec(NxD,NyD);
        Fv = [Fv zeros(size(Fv)); zeros(size(Fv)) Fv];
        Finv = [Finv zeros(size(Finv)); zeros(size(Finv)) Finv];
        Minv = Finv * MMinv * Fv;        
        M = [];
        if forwardOp == 1
            M = Finv * MM * Fv;
        end
    end

end
