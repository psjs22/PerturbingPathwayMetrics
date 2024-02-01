%%%%%%%% PERTURBED CONTRIBUTION METRICS. HYPOTHETICAL EXAMPLE 3.1 %%%%%%%%

clear all; close all; 

%%% COMPUTE MATRICES OF THE HYPOTHETICAL MODEL
mm = 1; % specify migratory model
[A,Ahat,n,c,s,mP,D,P,S,M,mD,mS,mM] = MigModel(mm);   

%%% COMPUTE UNPERTURBED CONTRIBUTION METRICS 
Phi = [1,2];    % full migratory routes
% compute all distict paths
PATH = DistinctPaths(n,s,Phi);  
%compute pathway metrics
[Cp,Ctilde] = PathwayMetrics(c,n,s,A,Ahat,Phi,PATH,mP);
% compute habitat metrics 
Ch = ones(1,4)*A(:,:,2)*A(:,:,1); 

%%% PERTURBATION STRUCTURE 
pD = sym(zeros(c,c,s)); % perturb demography 
syms d real; % act as delta 
Q=[0;1]; R=[0 1]; % structure of perturbation
% perturbation applies to D_{2,1}
pD(1:2,3:4,1) = D(1:2,3:4,1) + Q*d*R;
% all other D_{i,j} matrices are unperturbed
pD(1:2,1:2,1) = D(1:2,1:2,1); % D_{1,1}
pD(:,:,2) = D(:,:,2); % D_{1,2} and D_{2,2} 
% create perturbed versions of mD, A and Ahat 
pmD = sym(zeros(c*n,c*n,s));
pA = sym(zeros(c*n,c*n,s));
pAhat = sym(zeros(c*n,c*n,s));
for kk = 1:s
    % update demography in perturbed mD matrices
    for jj = 1:n
        E = zeros(n); E(jj,jj) = 1;     
        pmD(:,:,kk) = pmD(:,:,kk) + kron(E,pD(:,1+c*(jj-1):c*jj,kk));
    end 
    % update perturbed A and Ahat matrices
    pA(:,:,kk) = pmD(:,:,kk)*mM(:,:,kk); % movement happens first
    pAhat(:,:,kk) = pmD(:,:,kk)*mS(:,:,kk); % movement happens first 
end 

%%% COMPUTE PERTURBED CONTRIBUTION METRICS 
[pCP,pCPtilde] = PathwayMetrics(c,n,s,pA,pAhat,Phi,PATH,mP); 
% compute habitat contribution metrics 
pCh = ones(1,4)*pA(:,:,2)*pA(:,:,1);

%%% CHANGE IN CONTRIBUTION METRICS 
changeCP = pCP - Cp;
changeCPtilde = pCPtilde - Ctilde; 
changeCh = pCh - Ch; 



