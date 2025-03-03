
function plotEgTFMUQ(TFMUQ)

    magT = @(x) sqrt(mean(x.TTXUQ,3).^2 + mean(x.TTYUQ,3).^2);
    meanT = @(x) squeeze(mean(x,3));

    figure; 
    contourf(TFMUQ.X,TFMUQ.Y,magT(TFMUQ),'LineStyle','none'); 
    hold on;
    plots.quiverv2(TFMUQ.X,TFMUQ.Y,meanT(TFMUQ.TTXUQ),meanT(TFMUQ.TTYUQ),5e-3,...
            'Color','k','alpha',nan,'ratioHT',0.4,'linewidth',1);
    c=colorbar; 
    % set(c,'Ticks',[0 100 200]); 
    title(c,'Pa');
    plots.applyFormat;
    title('Traction magnitude');
end


