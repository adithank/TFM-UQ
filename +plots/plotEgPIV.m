
function plotEgPIV(out)

    calXY = 0.1628; % px -> um 
    magS = @(x) squeeze( sqrt( var(x.Usamp,0,3) + var(x.Vsamp,0,3) )); 

    [X,Y] = meshgrid(out.xvec,out.yvec);
    contourf(X,Y,magS(out)*calXY,'LineStyle','none'); % Total std. dev 
    hold on;
    
    plots.quiverv2(X,Y,out.UPost,out.VPost,5,'Color','w','alpha',nan,'ratioHT',0.4,'linewidth',1);
    c=colorbar; c.Ticks = [0 0.1 0.2];
    title(c,'\sigma_{PIVUQ} (\mum)');
    plots.applyFormat;
    title('PIV-UQ result');
end


