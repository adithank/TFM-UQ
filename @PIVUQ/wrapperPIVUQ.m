function output = wrapperPIVUQ(IMAGES, cfg_data, refFrames, sessFrames)


    % Assume refs are at the beginning
    Nrefs = length(refFrames);
    Ngrids = length(cfg_data.Deformation);
    Nsess = length(sessFrames);

    useGPU = cfg_data.Performance.UseGPU;
    
    % Number of bootstrap iterations 
    nB = [cfg_data.Deformation(:).Nreps];
    nB = min(nB(:));

    if length(unique(nB)) ~= 1
        warning('Multiple Nreps specified. Taking the minimum.');
        nB = min(nB(:));
    end

    % Send image to GPU
    if useGPU
        gpuDevice(1);
        IMAGES = gpuArray(IMAGES);
    end
    
    % Run without UQ once
    output = PIVUQ.wrapperPIVMP(IMAGES, cfg_data, refFrames, sessFrames);

    if nB == 1
        return;
    end

    % UQ run
    imR = IMAGES(:,:,refFrames);

    disp('Starting UQ run...');

    for i = 1 : Nsess

       if useGPU
           outUQ = UQGPU(imR,IMAGES(:,:,sessFrames(i)),cfg_data,nB);
           for igrid = 1 : Ngrids
                output(igrid).Usamp = outUQ(igrid).U; output(igrid).Vsamp = outUQ(igrid).V;
           end
       else
           outUQ = UQCPU(imR,IMAGES(:,:,sessFrames(i)),cfg_data,nB);

           for igrid = 1 : Ngrids
               output(igrid).Usamp = nan([size(output(igrid).U) nB]);
               output(igrid).Vsamp = nan([size(output(igrid).U) nB]);

                for ii = 1 : length(outUQ)
                    output(igrid).Usamp(:,:,ii) = outUQ{ii}(igrid).U;
                    output(igrid).Vsamp(:,:,ii) = outUQ{ii}(igrid).V;
               end
           end
           
       end
        
    end

    disp('UQ done');
    
    % Free-up GPU space
    if cfg_data.Performance.UseGPU
        gpuDevice(1);
    end
    

end



function outUQ = UQGPU(imR,imS,cfg_data,nB)

    N_RND = numel(imS) * nB;
    
    % Indices to drop in the image
    indS = randi(N_RND, [N_RND,1],'gpuArray');
    indS = setdiff(gpuArray.colon(1,N_RND),indS);
    
    imS = repmat(imS,[1 1 nB]);
    imS(indS) = 0;

    refFrames = 1;
    sessFrames = 1 + [1:nB];
    
    outUQ = PIVUQ.wrapperPIVMP(cat(3,imR,imS),cfg_data,refFrames,sessFrames);
        

end


function outUQ = UQCPU(imR,imSS,cfg_data,nB)

    % Only 2 images are supported with first as ref
    refFrames = 1; sessFrames = 2;
    Ngrids = length(cfg_data.Deformation);

    N_RND = numel(imSS);
    outUQ = cell(nB,1); % Init output cell arrays

    parfor ii = 1 : nB
    
        imS = imSS;
        rng('shuffle'); 
    
        % Indices to drop in the image
        indS = randi(N_RND, [N_RND,1]);
        indS = setdiff([1:N_RND],indS);
        imS(indS) = 0;
        
        output = PIVUQ.wrapperPIVMP(cat(3,imR,imS), cfg_data, refFrames, sessFrames);
        
        outUQ{ii} = output; 

    end


end

