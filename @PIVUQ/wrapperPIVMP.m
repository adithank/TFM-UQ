function output = wrapperPIVMP(IMAGES, cfg_data, refFrames, sessFrames, silentRun)

    assert(nargin>=4);

    if nargin == 4 
        silentRun = true; 
    end

    Ngrids = length(cfg_data.Deformation);
    output = struct();
     
    for igrid = 1:Ngrids

        if ~silentRun
            fprintf('Started pass %d ... \n', igrid);
        end 
        
        dum_cfg = cfg_data;
        dum_cfg.Deformation = dum_cfg.Deformation(igrid); % Per grid we modify a dummy config file that contains only the .Deformation needed for that pass
        
        
        if igrid == 1
            [X0,Y0,T0,U0,V0] = deal([]);
        else
            [X0,Y0,T0] = meshgrid(xvec,yvec,tvec); 
            U0 = U; V0 = V;
            % U0 = mean(U,4,'omitmissing'); V0 = mean(V,4,'omitmissing');
        end

        if ~silentRun
            PIVT = tic;
        end 
        
        % PIV step
        [xvec,yvec,tvec,U,V] = PIVUQ.PIVRic2024(IMAGES,dum_cfg,refFrames,sessFrames,X0,Y0,T0,U0,V0,[]);
        
        if ~silentRun
            PIVT = toc(PIVT); 
            fprintf('Finished pass %d in %1.2f s. \n', igrid, PIVT);
        end 
        
        output(igrid).xvec = xvec + dum_cfg.Deformation.wdw_size/2; % X,Y,T,R used to be 4-D matrices that took a lot of space in memory and provide little info. Now they are vectors. If you want the matrix form do [X,Y,T,R] =  meshgrid(xvec,yvec,tvec,repvec)
        output(igrid).yvec = yvec + dum_cfg.Deformation.wdw_size/2;
        output(igrid).U = U;
        output(igrid).V = V;
        
    end



end