function [verts] = quiverv2(X,Y,U,V,sc,options)
    % Adithan Kandasamy, 2023
    % Tested in Matlab R2023a
    % E.g. 
    % 1) quiverv2(X,Y,U,V)
    % 2) quiverv2(X,Y,U,V,sc=5,Color=[0 1 0],ratioHT=0.3, alpha=0.6, thetaH=20,lineWidth=1.5,useGPU=true);
    % 
    % INPUTS : 
    % X,Y   - Grid co-ordinates (matrix of size mxn)
    % U,V   - Vector values (matrix of size mxn)
    % sc    - Scaling of vectors (e.g. 0.1, 1 5 ,etc.)
    % Color - Character ('k','g',etc) or RGB vector ( [1 1 1], [0 1 0], etc.)
    % ratioTH - Ratio of head to tail (0 to 1, 0 for a line, 1 for no tail line)
    % alpha  - Transparency (0 to 1)
    % thetaH - Angle of head at the tip in degrees
    % linewidth - Width of head outline and tail line
    % useGPU - true / false, not much sppedup gain and GPU is slower the first time 

arguments 
    X double 
    Y double
    U double = X; % If nargin == 2
    V double = Y;
    sc double = 1 % Scaling of quivers
    options.Color (1,3) = [0 0 0] % RGB of both tail and head
    options.alpha (1,1) {mustBeNumeric} = 1 % Transparency [0,1]
    options.ratioHT (1,1) {mustBeNumeric} = 0.3 % Ratio of head to tail 
    options.thetaH (1,1) {mustBeNumeric} = 20 % Angle of head cone in degrees
    options.lineWidth (1,1) {mustBeNumeric} = 1.5 % Line width of the tail
    options.useGPU (1,1) = false % true for gpu use
end

if nargin <= 3 | isempty(X)
    X = [1:size(U,2)]; 
    Y = [1:size(U,1)];
    [X,Y] = meshgrid(X,Y);
end

C = options.Color;
fA = options.alpha;
rh = options.ratioHT;
thH = options.thetaH;
lw = options.lineWidth;

if ischar(C)
    C = C(1);
    C=bitget(find('krgybmcw'==C)-1,1:3);
end
% % Default parameters 
% rh = 0.3; % Ratio of head to total line length
% thH = 20; % Angle of glyph arrow head


% Check for size
assert(all(size(X) == size(Y)));
assert(all(size(X) == size(U)));
assert(all(size(Y) == size(V)));


% Scale the vectors
U = U*sc; V = V*sc;

% Get quiver positions
if options.useGPU
    X = gpuArray(X); Y = gpuArray(Y); U = gpuArray(U); V = gpuArray(V);
    % For GPU
    [p0x p0y pHx pHy pUHx pUHy pLHx pLHy] = arrayfun(@solveForQuiverPosition,X,Y,U,V,rh,thH);
else
    [p0x p0y pHx pHy pUHx pUHy pLHx pLHy] = arrayfun(@(X,Y,U,V) solveForQuiverPosition(X,Y,U,V,rh,thH),X,Y,U,V);
end

verts = struct();
verts.p0x = p0x; verts.p0y = p0y; 
verts.pHx = pHx; verts.pHy = pHy;
verts.pUHx = pUHx; verts.pUHy = pUHy;
verts.pLHx = pLHx; verts.pLHy = pLHy;

p0x = reshape(p0x,1,[]); p0y = reshape(p0y,1,[]);
pHx = reshape(pHx,1,[]); pHy = reshape(pHy,1,[]);
pUHx = reshape(pUHx,1,[]); pUHy = reshape(pUHy,1,[]);
pLHx = reshape(pLHx,1,[]); pLHy = reshape(pLHy,1,[]);

if isnan(fA)
    tt = sqrt(U(:).^2 +V(:).^2);
    thrs = prctile(tt,98,'all');
    tt(tt>=thrs) = max(tt(:));
    fA = tt./max(tt(:));
    patch([p0x;(pUHx+pLHx)/2],[p0y;(pUHy+pLHy)/2],C,'faceAlpha','flat','edgeColor',C,'edgeAlpha','flat','lineWidth',lw,...
                            'faceVertexAlphaData',fA);
    patch([pUHx;pHx;pLHx],[pUHy;pHy;pLHy],C,'EdgeColor',C,'FaceAlpha','flat','EdgeAlpha','flat','lineWidth',lw,...
                            'faceVertexAlphaData',fA);
else 
    patch([p0x;(pUHx+pLHx)/2],[p0y;(pUHy+pLHy)/2],C,'faceColor',C,'faceAlpha',fA,'edgeColor',C,'edgeAlpha',fA,'lineWidth',lw);
    patch([pUHx;pHx;pLHx],[pUHy;pHy;pLHy],C,'faceColor',C,'EdgeColor',C,'FaceAlpha',fA,'EdgeAlpha',fA,'lineWidth',lw);
end

% Plot lines 
% clf;
% line([p0x;pHx],[p0y;pHy],'Color',[C fA],'lineWidth',lw);

% line([pUHx;pHx],[pUHy;pHy],'Color','k','lineWidth',1.5);
% line([pLHx;pHx],[pLHy;pHy],'Color','k','lineWidth',1.5);
% hold on; 




% hold off;
% patch([pUHx;pHx;pLHx],[pUHy;pHy;pLHy],C,'EdgeColor',C,'FaceAlpha',fA,'EdgeAlpha','interp','lineWidth',1,...
%     'faceVertexAlphaData',normalize(sqrt(U(:).^2 +V(:).^2)));



% clf;
% lirne(lineX1,lineY1,'Color','k','lineWidth',2);
% line(lineX2,lineY2);
% line(lineX3,lineY3);
% patch([O(1)+lb O(1)+l O(1)+lb],[O(2)+h O(2) O(2)-h],'k','lineWidth',2);
% axis([-1 1 -1 1])

end


function [p0x p0y pHx pHy pUHx pUHy pLHx pLHy] = solveForQuiverPosition(X,Y,U,V,rh,thH)
    % Solve for body positions 
    %       pUH
    % p0 ----> pH
    %       pLH
    %   lb - Length of body 
    %   lh - Length of head 
    %   h  - height of arrow head 


    if isnan(X) | isnan(Y) | isnan(U) | isnan(V)
        p0x = nan; p0y = nan; pHx = nan; pHy = nan; pUHx = nan; pUHy = nan; pLHx = nan; pLHy = nan;
        return;
    end

    l = sqrt(U.^2 + V.^2);
    alp = atan2(V,U);

    lh = l*rh;
    lb = l - lh; 
    h = lh*tand(thH);

    p0x = X; p0y = Y;
    % pH = [l 0];
    % pUH = [lb h];
    % pLH = [lb -h];
    % 
    % RT = [cos(alp) -sin(alp); sin(alp) cos(alp)];
    % pH = RT*pH';   pHx=pH(1)+X;   pHy=pH(2)+Y;
    % pUH = RT*pUH'; pUHx=pUH(1)+X; pUHy=pUH(2)+Y;
    % pLH = RT*pLH'; pLHx=pLH(1)+X; pLHy = pLH(2)+Y;

    % Changed for GPU 
    pHxT = l; pHyT = 0;
    pUHxT = lb; pUHyT = h;
    pLHxT = lb; pLHyT = -h;

    RT11 = cos(alp); RT12 = -sin(alp); RT21 = sin(alp); RT22 = cos(alp);

    pHx=RT11*pHxT+RT12*pHyT+X;   pHy=RT21*pHxT+RT22*pHyT+Y;
    pUHx=RT11*pUHxT+RT12*pUHyT+X;   pUHy=RT21*pUHxT+RT22*pUHyT+Y;
    pLHx=RT11*pLHxT+RT12*pLHyT+X;   pLHy=RT21*pLHxT+RT22*pLHyT+Y;

end
