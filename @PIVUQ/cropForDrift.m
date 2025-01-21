

function [refIm, sessIm] = cropForDrift(refIm, sessIm, udis, vdis)

    % Drift is calculated as drift of sess image from ref image
    % (i.e.) ref image is cropped in place 
    % sess image is drift corrected and cropped
    % [udis, vdis] is 2D Xand Y drift component
    % 
    % ref or sess can be empty vector 

    % Only corrects upto integer
    udis = round(udis); vdis = round(vdis);
    
    % X
    if udis > 0
        sessx = [udis+1:size(sessIm,2)];
        refx = [1:size(refIm,2)-udis];
    elseif udis < 0
        sessx = [1:size(sessIm,2)+udis];
        refx = [1-udis:size(refIm,2)];
    elseif udis == 0
        sessx = [1:size(sessIm,2)];
        refx = [1:size(refIm,2)];
    end

    % Y
    if vdis > 0
        sessy = [vdis+1:size(sessIm,1)];
        refy = [1:size(refIm,1)-vdis];
    elseif vdis < 0
        sessy = [1:size(sessIm,1)+vdis];
        refy = [1-vdis:size(refIm,1)];
    elseif vdis == 0
        sessy = [1:size(sessIm,1)];
        refy = [1:size(refIm,1)];
    end

    refIm = refIm(refy,refx);
    sessIm = sessIm(sessy,sessx);


end