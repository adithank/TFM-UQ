% Using kronecker product (outer product) to keep things simple
function [Fvec,Finv] = genDFTmatricesFor2DVec(m,n)
    Fvec = kron(dftmtx(m),dftmtx(n));
%     Finv = inv(Fvec);    
    Finv = Fvec \ speye(size(Fvec));
end







%% Old pedantic method
% % Copied from https://www.mathworks.com/matlabcentral/answers/153214-how-to-create-2d-dft-matrix-to-transform-a-vectorized-2d-image
% 
% function [Fvec,Finv] = genDFTmatricesFor2DVec(m,n)
%     Fvec = zeros(m*n,m*n);
%     for i=1:m*n
%         inpvec = zeros(m*n,1);
%         inpvec(i) = 1;
%         inpvec = reshape(inpvec,m,n);
%         %finding out one col of the transformation. since i/p and o/p basis are
%         %canonical, nothing further is required
%         opvec = fft2(inpvec); 
%         opvec = reshape(opvec,[],1);
%         Fvec(:,i) = opvec;
%     end
% Finv = inv(Fvec);    
% end

% Inverse of Fvec is IDFT matrix


% Tried this, didnot work. Shape looks correct, but it doesnt seem to give
% correct answers. Getting concentration at (1,1).
% THINK PROBLEM IS THAT IM NOT FFTSHIFTING
% F1dn = dftmtx(20);
% for i = 1 : 20
% FLTmp(i,:) = [zeros(1,(i-1)*20) F1dn(i,:) zeros(1,(20-i)*20)];
% end
% FLeft = repmat(FLTmp,20,1);
% FRight = zeros(20*20,20*20);
% for k = 1 : 20
% clear FRTmp;
% for i = 1 : 20
% for j = 1 : 20
% FRTmp(i,((j-1)*20)+i) = F1dn(j,k);
% end
% end
% FRight((k-1)*20+1:(k-1)*20+20,:) = FRTmp;
% end
% 
% Fvec = FRight * FLeft; 