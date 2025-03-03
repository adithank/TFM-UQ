classdef PIVUQ < handle


properties (Access = public)
    cfg_data = []; 
    IMAGES = [];
    refFrames = [];
    sessFrames = [];
    drift = [];

    output = []; % Results
end



methods (Access = public)

    function obj = correctForDrift(obj)

        assert(isscalar(obj.refFrames),'Only one refFrame is allowed');
        assert(~isempty(obj.refFrames) | ~isempty(obj.sessFrames),'No ref and sessframes specified');

        dum_cfg = obj.loadJsonConfig('defaultPIVconfig.json');
        dum_cfg.Deformation.Nreps = 1;
        dum_cfg.Deformation.wdw_size = fix(size(obj.IMAGES,1));
        dum_cfg.Deformation.wdw_spacing = fix(size(obj.IMAGES,1)/2);

        [~,~,~,udis,vdis] = obj.PIVRic2024(obj.IMAGES,dum_cfg,obj.refFrames,obj.sessFrames);
        
        [refIm] = obj.cropForDrift(obj.IMAGES(:,:,obj.refFrames), [], udis, vdis);
        sessIm = nan([size(refIm), size(obj.IMAGES,3)-1]);
        
        for ii = 1 : size(sessIm,3)
            [~,tt] = obj.cropForDrift([], obj.IMAGES(:,:,obj.sessFrames(ii)), udis, vdis);
            sessIm(:,:,ii) = tt;
        end
        obj.IMAGES = cat(3,refIm,sessIm);
        obj.drift = struct();
        obj.drift.udis = udis; obj.drift.vdis = vdis; 
    end

end

methods (Access = public)

   % Constructor
   function obj = PIVUQ(configPath)

       if nargin == 1
           obj.cfg_data = obj.loadJsonConfig(configPath);
       else
           obj.cfg_data = obj.loadJsonConfig('defaultPIVConfig.json');
       end
   end
       
   function obj = runPIVMPUQ(obj)

        obj.output = obj.wrapperPIVUQ(obj.IMAGES, obj.cfg_data, obj.refFrames, obj.sessFrames);

   end

   function obj = postProcess(obj)
       obj.output = obj.postProcessPIVUQ(obj.output);
       
       % Remove drift by median subtraction 
       for ii = 1: length(obj.output)
            obj.output(ii).UPost = obj.output(ii).UPost - median(obj.output(ii).UPost(:),'omitmissing');
            obj.output(ii).VPost = obj.output(ii).VPost - median(obj.output(ii).VPost(:),'omitmissing');
       end
   end

end



methods (Static)
    
    [cfg_data] = loadJsonConfig(configfile);
    
    output = wrapperPIVUQ(IMAGES, cfg_data, refFrames, sessFrames, silentRun);

    output = wrapperPIVMP(IMAGES, cfg_data, refFrames, sessFramesm, silentRun);

    [xvec,yvec,tvec,U,V] = PIVRic2024(IMAGES,dum_cfg,refFrames,sessFrames,X0,Y0,T0,U0,V0,C0);
    
    [refIm, sessIm] = cropForDrift(refIm, sessIm, udis, vdis);

    output = postProcessPIVUQ(output);

end


end