function quiverv2_UQcone(X,Y,Usamp,Vsamp,nSD,sc,options)

arguments 
    X
    Y
    Usamp
    Vsamp
    nSD double = 1 % Scaling of standard deviation of mag and dir.
    sc double = 1 % Scaling of quivers
    options.EdgeColor = 'k';
    options.FaceColor = 'k';
    options.EdgeAlpha (1,1) {mustBeNumeric} = 0.8 % Transparency [0,1]
    options.FaceAlpha (1,1) {mustBeNumeric} = 0.3 % Transparency [0,1]
    options.ratioHT (1,1) {mustBeNumeric} = 0.3 % Ratio of head to tail 
    options.lineWidth (1,1) {mustBeNumeric} = 1.5 % Line width of the tail
    options.angularOnly (1,1) = false; % true to plot equal length gylphs
    options.plotMean (1,1) = false; % Plot mean gylphs
    options.useGPU (1,1) = false % true for gpu use
end

fC = options.FaceColor;
eC = options.EdgeColor;
eA = options.EdgeAlpha;
fA = options.FaceAlpha;
rh = options.ratioHT;
lw = options.lineWidth;
angOnly = options.angularOnly;

% Re-scale variables
% U = U *sc; V = V * sc;
if ndims(Usamp) > 3 
    Usamp = squeeze(Usamp); Vsamp = squeeze(Vsamp);
end

Usamp = Usamp * sc; Vsamp = Vsamp * sc;

% Plot mean gylph

if options.plotMean & ~angOnly 
    U = trimmean(Usamp,10, 3); V = trimmean(Vsamp,10,3);
    [verts] = quiverv2(X,Y,U,V); 
end


% Convert components to dtheta and dm
[m,dm,th,dth] = convertcompToAngleMag(Usamp,Vsamp,nSD);

if angOnly
    m = ones(size(m))*sc;
end

% Delete big fans
dth(abs(dth)>3*pi/4) = nan;

% Get cone positions
% p0x = verts.p0x; p0y = verts.p0y;
% alp = atan2(V,U);


if isgpuarray(X)
    return;
else
    [p0X,p0Y,p1X,p1Y,p3X,p3Y,p5X,p5Y, ...
            p1HX,p1HY,p3HX,p3HY,p5HX,p5HY, ...
                p1LX,p1LY,p3LX,p3LY,p5LX,p5LY] = arrayfun( @(X,Y,m,dm,th,dth) solveForConePosition(X,Y,m,dm,th,dth,rh),...
                                            X,Y,m,dm,th,dth);

    % [p0X,p0Y,p1X,p1Y,p3X,p3Y,p5X,p5Y] =  solveForConePosition(X(1),Y(1),m(1),dm(1),th(1),dth(1),ratioTH);

    % yLC = xLC; yHC = xLC; xHC = xLC;
    % for ii = 1 : numel(X)
    %     [xLCtt yLCtt xHCtt yHCtt] = solveForConePosition(X(ii),Y(ii),U(ii),V(ii),m(ii),dm(ii),th(ii),dth(ii),alp(ii),p0x(ii),p0y(ii),nCone);
    %     xLC(:,ii) = xLCtt; yLC(:,ii) = yLCtt; xHC(:,ii) = xHCtt; yHC(:,ii) = yHCtt;
    % end
end

hold on; 
patch([p0X(:)';p1X(:)';p3X(:)';p5X(:)';p0X(:)'],[p0Y(:)';p1Y(:)';p3Y(:)';p5Y(:)';p0Y(:)'],fC,'edgeColor',eC,'faceColor',fC,'faceAlpha',fA,...
                        'edgeAlpha',eA,'linewidth',lw);

if angOnly
    patch([p0X(:)'; p3X(:)'],[p0Y(:)'; p3Y(:)'],eC,'edgeColor',eC,'edgeAlpha',eA,'lineWidth',lw)

else
    filler = nan(size(p1HX));
    
    % Uncertainty in magnitude
    patch([p1HX(:)';p3HX(:)';p5HX(:)';filler(:)'],[p1HY(:)';p3HY(:)';p5HY(:)';filler(:)'],fC,'edgeColor',eC,'faceColor',fC,'faceAlpha',0,...
                            'edgeAlpha',eA,'linewidth',lw);
    patch([p1LX(:)';p3LX(:)';p5LX(:)';filler(:)'],[p1LY(:)';p3LY(:)';p5LY(:)';filler(:)'],fC,'edgeColor',eC,'faceColor',fC,'faceAlpha',0,...
                            'edgeAlpha',eA,'linewidth',lw);

end


% hold on; quiver(X,Y,U,V,0,'w','LineWidth',1.5)


% clf;
% lirne(lineX1,lineY1,'Color','k','lineWidth',2);
% line(lineX2,lineY2);
% line(lineX3,lineY3);
% patch([O(1)+lb O(1)+l O(1)+lb],[O(2)+h O(2) O(2)-h],'k','lineWidth',2);
% axis([-1 1 -1 1])

end

function [m,dm,th,dth] = convertcompToAngleMag(Usamp,Vsamp,nSD)
    magSamp = sqrt(Usamp.^2 + Vsamp.^2); 
    % m = trimmean(magSamp,10,3);
    Um = trimmean(Usamp,10,3); Vm = trimmean(Vsamp,10,3);
    m = sqrt( Um.^2 + Vm.^2 );
    dm = squeeze( sqrt( trimmean( (magSamp - m).^2 ,10,3) ) ); % Trimmed std dev
    dm = dm * nSD;
    
    % Angle of mean w.r.t x-axis
    % phi = atan2()

    thsamp = atan2(Vsamp,Usamp);
    % th = circ_mean(thsamp,[],3);
    th = atan2(Vm,Um);
    [~,dth] = circ_std(thsamp - th,[],[],3);
    dth = dth * nSD;
end

function [p0X,p0Y,p1X,p1Y,p3X,p3Y,p5X,p5Y, ...
            p1HX,p1HY,p3HX,p3HY,p5HX,p5HY, ...
                p1LX,p1LY,p3LX,p3LY,p5LX,p5LY] = solveForConePosition(X,Y,m,dm,th,dth,ratioTH)
    % Solving for 0 1 3 5 0 points in the paper, scaled to magnitude not
    % area (i.e. l = m)
    % p_i denote vertices

    % Solve for la, lb, h
    % lb = sqrt( m./( tan( dth/2 ) .* (ratioTH+1) )  ); 
    % la = ratioTH .* lb; 
    % h = m/(la + lb); 

    % Set m = la + lb
    [la,lb,h] = solveForConeParams(m,ratioTH,dth);
    
    % Solve for wingDx wingDy
    % wingDx = ratioTH.*la; wingDy = ratioTH.*h; 

    p0X = X;    poY = Y;
    % Base cone with angular uncertainty
    p1X = lb;   p1Y = -h; 
    p3X = lb+la;p3Y = 0;
    p5X = lb;   p5Y = h;

    % Uncertainty in magnitude
    [la,lb,h] = solveForConeParams(m+dm/2,ratioTH,dth);
    p1HX = lb;   p1HY = -h; 
    p3HX = lb+la;p3HY = 0;
    p5HX = lb;   p5HY = h;

    [la,lb,h] = solveForConeParams(m-dm/2,ratioTH,dth);
    p1LX = lb;   p1LY = -h; 
    p3LX = lb+la;p3LY = 0;
    p5LX = lb;   p5LY = h;

    % Rotation transformation
    RT11 = cos(th); RT12 = -sin(th); RT21 = sin(th); RT22 = cos(th);
    
    p0X = X; p0Y = Y;
    % [p0X,p0Y] = rotatePoints(p0X,p0Y);
    [p1X,p1Y] = rotatePoints(p1X,p1Y,RT11,RT12,RT21,RT22);
    [p3X,p3Y] = rotatePoints(p3X,p3Y,RT11,RT12,RT21,RT22);
    [p5X,p5Y] = rotatePoints(p5X,p5Y,RT11,RT12,RT21,RT22);

    [p1HX,p1HY] = rotatePoints(p1HX,p1HY,RT11,RT12,RT21,RT22);
    [p3HX,p3HY] = rotatePoints(p3HX,p3HY,RT11,RT12,RT21,RT22);
    [p5HX,p5HY] = rotatePoints(p5HX,p5HY,RT11,RT12,RT21,RT22);

    [p1LX,p1LY] = rotatePoints(p1LX,p1LY,RT11,RT12,RT21,RT22);
    [p3LX,p3LY] = rotatePoints(p3LX,p3LY,RT11,RT12,RT21,RT22);
    [p5LX,p5LY] = rotatePoints(p5LX,p5LY,RT11,RT12,RT21,RT22);
    
    p1X = p1X + p0X; p1Y = p1Y + p0Y;
    p3X = p3X + p0X; p3Y = p3Y + p0Y;
    p5X = p5X + p0X; p5Y = p5Y + p0Y;

    p1LX = p1LX + p0X; p1LY = p1LY + p0Y;
    p3LX = p3LX + p0X; p3LY = p3LY + p0Y;
    p5LX = p5LX + p0X; p5LY = p5LY + p0Y;

    p1HX = p1HX + p0X; p1HY = p1HY + p0Y;
    p3HX = p3HX + p0X; p3HY = p3HY + p0Y;
    p5HX = p5HX + p0X; p5HY = p5HY + p0Y;
end

function [la,lb,h] = solveForConeParams(m,ratioTH,dth)
    la = ratioTH.*m;
    lb = m - la;
    h = lb.*tan(dth/2);
end

function [pX,pY] = rotatePoints(pXT,pYT,RT11,RT12,RT21,RT22)
        pX = RT11*pXT+RT12*pYT;      pY= RT21*pXT+RT22*pYT;
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



% function [m,dm,th,dth] = convertcompToAngleMag(U,V,Usamp,Vsamp,nSD)
%     magSamp = sqrt(Usamp.^2 + Vsamp.^2); 
%     m = trimmean(magSamp,10,3);
%     dm = squeeze( sqrt( trimmean( (magSamp - m).^2 ,10,3) ) ); % Trimmed std dev
%     dm = dm * nSD;
% 
%     thsamp = atan2(Vsamp,Usamp);
%     th = circ_mean(thsamp,[],3);
%     [~,dth] = circ_std(thsamp - th,[],[],3);
%     dth = dth * nSD;
% end