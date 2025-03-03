function [alpV,betV,TT,it] = MCMCGibbs_TFMUQ(M,UU,alp0,bet0,L,SigPIV,params,prevIts,saveIts)

% Inputs :
%
% M          - Elastostatic operator M : T -> U 
% [UU]        - Displacement field values in vector form [-U-; -V- ; -W-]; 
% alp0       - Intial value for alpha (Prior) hyperparameter
% bet0       - Initial value for beta (Global error) hyperparameter( inf to
%               skip beta smapling )
%
%   IF not provided, MLE estimate is used to initialize
%
% L          - Precision matrix of traction Prior distribution
% SigPIV     - Variance matrix of U from PIV-UQ 
% params     - Structure with fields : maxIt (maximum Gibbs iterations),
%               burnIt (Burn-in iterations), plotT (true to plot each
%               iteration result)
% prevIts    - Structure with fields : it (iteration count), alpV (Array of
%               sampled alp values), betV (Vector of sampled beta values),
%               TT (Vectpr of sampled posterior traction vectors "t")  
%
% saveIts    - Structure with fields : saveFreq (Iteration freq to save, 0 for no save), 
%               savePath (Path with filename)
%
% OUTPUT : 
% alpV       - Vector of sampled alpha hyper-parameter after burnin
% betV       - Vector of sampled beta  hyper-parameter after burnin
% TT         - Vector of sampled traction vectors
% it         - Last iteration of Gibbs sampler


% Parameters for the sampler
maxIt  = params.maxIt;
burnIt = params.burnIt;
plotT = params.plotT; 

% Resume iterations if not empty
if ~isempty(prevIts)
    it     = prevIts.it;
    alpV  = prevIts.alpA; alp = alpV(end);
    betV  = prevIts.betA; bet = betV(end);
    TT   = prevIts.TT; 
else
    it = 1;
    alpV = []; betV = []; TT = []; 
    alp = alp0; bet = bet0;
end


% Fixed parameters for hyper-priors (Gamma distribution)
thB = 1; phiB = 1e-5; 
thA = 1; phiA = 1e-5;

n = numel(UU);

U = speye(size(M));


statCounter = 0; % Counter for stationary chain, allow upto 5 retries



% try
    while it < (maxIt+burnIt)
        
        %%% %%% %%% %%% %%% NICE status printing here %%% %%% %%% %%%

        % Likelihood precision matrix, Omega^-1 = bet^-1 I + SigPIV , I is identity
        % OmeInv = inv((1/bet)*U + SigPIV);  %%% Use ldivide instead of inv 
        OmeInv = ((1/bet)*U + SigPIV) \ speye(size(SigPIV,2));
        
        % Posterior covariance
        cov = real( (  M.'* OmeInv *M + alp * L)) \ speye(size(M,2)) ;
        cov = 0.5*(cov + cov.'); % Symmetrisize
        % cov(cov<0) = 0; % cov 

        tMean =  real( cov * M.'* OmeInv *UU ); % Conditional posterior mean
        
        try
            tt = mvnrnd(tMean,cov,1); % Sample a traction posterior 
        catch e
            keyboard; 
            % IF unable to sample
            statCounter = statCounter + 1;
            if statCounter >= 5 
                warning(' MCMC chain is unable to find the next iteration. Terminating ... '); 
                return;
            else        
                % If fails, try to sample again  
                continue; 
            end
        end

        tt = real(tt.'); 

        alp = gamrnd( n/2 + thA, 1./( 0.5 * (tt.' * L * tt) + phiA) );

        if ~isinf(bet) 
            
            betStatCounter = 0;

            % Direct sampling
            [betT] =  sampleBetawUQ(UU-real(M*tt),diag(SigPIV),phiB,thB);
            
            if betT == - 1 % Error in direct sampling
               betStatCounter = betStatCounter + 1; 

               if betStatCounter > 3 
                   return;
               end

               continue;   
            end
            
            bet = betT; 
        end
        
       it=it+1;

       % If sampling is good and store values after burn-in 
       if it > burnIt
           ind = it - burnIt;

           TT(:,ind) = tt;
           alpV(ind) = alp;
           betV(ind) = bet;
            
           fprintf(' It - %d; Sqrt(alp) - %2.2f sqrt(bet) - %2.2f \n',ind, sqrt(1./alp),sqrt(1./bet));
           
           if ~isempty(saveIts)
                
                if mod(it,saveIts.freq) == 0
                    save(saveIts.savePath, 'alpV','betV','it','TT');
                end

           end
            
           if plotT 
                figure(101);
                subplot(3,1,1)
                hold on;
                plot(it,sqrt(mean( ( real(UU - M * tt) ).^2 ,'omitmissing') ) , 'b-*'); % RMSE
                ylabel('Residual RMSE ');
                subplot(3,1,2)
                hold on;
                semilogy(it, sqrt(1./bet), 'k-o');
                ylabel('sqrt(1/\beta)');
                subplot(3,1,3)
                hold on;
                semilogy(it, sqrt(1./alp), 'r-o');
                ylabel('Prior variance, (1/\alpha)^2');
                xlabel('Iteration');
                drawnow;
    
                % figure(102)
                % clf;
                % imagesc(real( reshape(tt(1:end/2),size(XT)) ));
                % title('T_x realization');
                % colorbar;
                % drawnow;
            end
        end

    end
    
    
end



function [bet] =  sampleBetawUQ(res,Zii,phiB,thB)

% INPUTS 
% 
% res - Residual u - M t
% Zii - Diagonal vector of SigPIV
% phi, th - Hyper-prior constant parameters 

    % log p(beta | uh ... )

    logyyfun = @(bet) -( -0.5*(res(:).'*diag(1./(1./bet + Zii))*res(:)) - phiB*bet ...
                            + 0.5*sum(log( 1./( 1./bet + Zii )))...
                                +(thB-1)*log(bet) );
                    
	options = optimoptions('fminunc','Display','none','Algorithm','quasi-newton');
    [betHat,~,exitflag,~,~,hessian] = fminunc(logyyfun,1,options); % betHat is MLE estimate 

    if (exitflag~=1&exitflag~=2) % Check for convergence
        disp(exitflag);
        bet = -1; return;
    end
    
    % Construct a strictly positive discrete domain for beta around MLE
    % estimate (+/-3 std) with 100 points

    betArr = linspace(betHat-sqrt(1./hessian)*3,betHat+sqrt(1./hessian)*3,100);
    betArr(betArr<=0)=[];

    Logyy = nan(length(betArr),1);
    for i = 1 : length(betArr)
        Logyy(i) = -logyyfun(betArr(i));
    end

    % pdf
    pp = (exp((Logyy) - max(Logyy)));
    pp = pp./sum(pp);

    % cdf and inverse sampling
    ii=[];for j = 1 : 10
        [~,ii] = min(abs(rand-cumsum(pp)));
    end

    bet = betArr(ii);

end