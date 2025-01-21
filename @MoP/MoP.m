classdef MoP

    methods(Static = true)
    
        
            [M,Minv] = TfmMatrices2DFiniteThickLinButler(XD,YD,XT,YT,E,nu,h,realDom,forwardOp,Zero0thmode);
        
            [Fvec,Finv] = genDFTmatricesFor2DVec(m,n)


        end
    
end



