function [Cp,Cptilde,AAhat,AAhatP] = ... 
    PathwayMetrics(c,n,s,A,Ahat,Phi,PATH,mP)
% Compute pathway metrics for PATH
% Cp contains the subpopulation pathway contribution metrics. See (2.12) in
% the paper https://doi.org/10.1016/j.ecolmodel.2022.110056
% Cptilde contains the metapopulation pathway contribution metrics. See
% (2.13) in the paper. 
% AAhat contains the product of the seasonal and/or seasonal survival
% matrices for the specified path. See (2.11) in the paper. 
% AAhatP contains the product of the seasonal and/or seasonal survival
% matrices along with the matrices containing the proportion of the 
% population using the specified path. See (2.13) and next two equations 
% in the paper.

AAhat = zeros(c*n,c*n,size(PATH,1)); % to store product of A and Ahat matrices for each distinct path
AAhatP = zeros(c*n,c*n,size(PATH,1)); % to store product of A, Ahat and \bP matrices for each distinct path
Cp = zeros(size(PATH,1),c*n); % to store subpopulation pathway metrics for each distinct path
Cptilde = zeros(size(PATH,1),c*n); % to store metapopulation pathway metrics for each distinct path

for cc = 1:size(PATH,1) % work through all distinct paths
    
% CONSTRUCT THE REQUIRED PRODUCT OF A AND Ahat MATRICES FOR EACH PATH
    AAhat(:,:,cc) = eye(c*n,c*n);
    AAhatP(:,:,cc) = eye(c*n,c*n); 
    for kk = 1:s% anniversary season labelled as 1
        % If the path taken in season kk IS specified
        if ismember(kk,Phi) %PATH(cc,kk)~=0 && PATH(cc,kk+1)~=0
            E = zeros(n); E(PATH(cc,kk+1),PATH(cc,kk)) =1;
            AAhat(:,:,cc) = (Ahat(:,:,kk).*kron(E,ones(c,c)))...
                *AAhat(:,:,cc);
            AAhatP(:,:,cc) = (Ahat(:,:,kk).*kron(E,ones(c,c)))...
                *kron(eye(n),mP((PATH(cc,kk+1)-1)*c+1:PATH(cc,kk+1)*c,...
                        (PATH(cc,kk)-1)*c+1:PATH(cc,kk)*c,kk))...
                *AAhatP(:,:,cc);
        % If the path taken in season kk is NOT specified
        else % path is unspecified in season kk
            AAhat(:,:,cc) = A(:,:,kk)*AAhat(:,:,cc);
            AAhatP(:,:,cc) = A(:,:,kk)*AAhatP(:,:,cc);
        end 
    end
    
% CALCULATE PATHWAY METRICS
    Cp(cc,:) = ones(1,c*n)*AAhat(:,:,cc);
    Cptilde(cc,:) = ones(1,c*n)*AAhatP(:,:,cc);
end 
