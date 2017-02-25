function X = myLinCG(R)
%% myLinCG solves the linear equation M * X = R for X
% Inputs:
%   R: known residuals
% Outputs:
%   X: solution

global gLG gLH gZ gMu gMaxitIn

% initialization
n = size(gLG,1);
nR = size(R,2);
X = sparse(n,nR);
flag = zeros(nR,1);                                 

tol=1e-3;
maxIt = gMaxitIn;

M = diag(gLG)+ gMu*diag(gLH) + gZ/sqrt(n);
M = sparse(M);

for i=1:nR;
    x=X(:,i);
    b=R(:,i);
    bnrm2 = norm( b );
    if  ( bnrm2 == 0.0 ), 
        X(:,i)=zeros(size(x));
        flag(i)=1;
        continue
    end
    regularAx = gLG*x + gMu*(gLH*x) + sum(x)*gZ/sqrt(n);
    r = b - regularAx;
    error = norm( r ) / bnrm2;
    if ( error <= tol ) 
        flag(i)=1;
        continue
    end
    % begin iteration
    for iter = 1:maxIt                        
        z  = r./M;
        rho = r'*z;
        % direction vector
        if iter > 1                   
            beta = rho / rho_1;
            p = z + beta*p;
        else
            p = z;
        end
		regularAp = gLG*p + gMu*(gLH*p) + sum(p)*gZ/sqrt(n);
        alpha = rho / (p'*regularAp );
        % update approximation vector
        x = x + alpha * p;                    
        % compute residual
        r = r - alpha * regularAp;                      
        error = norm( r ) / bnrm2;
        % check convergence
        if ( error <= tol )
            flag(i) = 1;
            disp('inner iteration converges');
            break
        end 
        rho_1 = rho;
    end
    X(:,i)=x;
end
end
