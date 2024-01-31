%%%%%%%%%%%%%%%%% PERTURBING PATHWAY CONTRIBUTION METRICS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% MONARCH BUTTERFLY EXAMPLE %%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
% Compute matrices for model mm in {1,2,3,4,5,6,7,8,9} 
mm = 6; % specify migratory model
[A,Ahat,n,c,s,mP,D,P,S,M,mD,mS,mM] = MigModel(mm);

% Compute the unperturbed PCMs for each individual path in season 6
Phi = 6; % SPECIFYING A PATH IN JUST SEASON 6
PATH = zeros(4,s+1); %4 distict paths associated with season 6
% populate PATH1 
PATH(1:2,6)=2; PATH(3:4,6)=3; % all valid start points
PATH([1,3],7)=3; PATH([2,4],7)=4; % all valid end points 

Cp = zeros(size(PATH,1),c*n);
Cptilde = zeros(size(PATH,1),c*n);
ChPtilde = zeros(size(PATH,1),c*n); 
% compute all unperturbed pathway contribution metrics
[Cp(:,:),Cptilde(:,:),Cptildej,Cptildei] = ...
        PathwayMetrics(c,n,s,A,Ahat,Phi,PATH(:,:),mP);    
% compute habitat metrics from Cptilde
ChPtilde(1,1:c*n) = ones(1,size(PATH,1))*Cptilde(:,:);


%% structural matrices and constants in the sensitivity formula

% independent of path
K = zeros(c*n); % vec permutation matrix K_{c,n}
for i=1:c
    for j=1:n
        E = zeros(c,n);
        E(i,j) = 1; 
        K = K + kron(E,E');
    end
end

X1 =  kron(kron(eye(n),K),eye(c)); % constant X_1
X4 = kron(eye(c*n),ones(1,c*n)); % constant X_4

X2k = zeros((c*n)^2,(c*n)^2,s); % constant X_{2,k}
for kk=1:s
    X2k(:,:,kk) = kron(eye(c*n),mD(:,:,kk));
end

Ec44 = zeros(c); Ec44(4,4)=1; %matrix E_{c,ii} i=4
vecEc44 = reshape(Ec44,c*c,1);

Ec55 = zeros(c); Ec55(5,5)=1; %matrix E_{c,ii} i=5
vecEc55 = reshape(Ec55,c*c,1);

% constant for a specified path in season 6
Z1k = zeros((c*n)^2,(c*n)^2,s); %constant Z_{1,k}
for kk = 1:s
    if kk==6    % season 6 has a specified path
        Z1k(:,:,kk) = kron(mS(:,:,kk)',eye(c*n));
    else        % all other seasons do not have a specified path
        Z1k(:,:,kk) = kron((mS(:,:,kk).*mP(:,:,kk)),eye(c*n));
    end
end

Z46 = eye(c*n*c*n); %constant Z_{4,6} when there is a specified path in k=6

Y51 = eye(c*n); % to store constant Y^5_1
Y127 = Y51; % to store constant Y^12_7
for kk=1:5   
    Y51 = A(:,:,kk)*Y51;
    Y127 = A(:,:,kk+6)*Y127;    %only goes up to season 11
end
Y127 = A(:,:,12)*Y127; 
kronYY = kron(Y51',Y127); % (Y^5_1)^T \otimes Y^12_7 

% \Upsilon^5_1 = Y^5_1 as none of these seasons have a specified path and
% \Upsilon^12_7 = Y^12_7 for same reasons
U51 = Y51;
U127 = Y127; 
kronUU = kronYY; 

vecmP6 = reshape(mP(:,:,6),c*n*c*n,1);
vecmS6 = reshape(mS(:,:,6),c*n*c*n,1);

vecIn = reshape(eye(n),n*n,1);

% Diagonal matrices used in the calculation of elasticities. 
diagCp =  zeros(20,20,size(PATH,1)); % a matrix associated with each path
diagCptilde = diagCp;
diagChabitat = diagCp;
for ii = 1:size(PATH,1)
    diagCp(:,:,ii) = diag(Cp(ii,:));
    diagCptilde(:,:,ii) = diag(Cptilde(ii,:));
    diagChabitat(:,:,ii) = diag(ChPtilde(ii,:));
end 

%% SENSITIVITY OF CONTRIBUTION METRICS WITH RESPECT TO DEMOGRAPHIC RATES
% MIGRATION HAPPENS BEFORE DEMOGRAPHY => PATH(i,7) indicates which
% habitat's life rates will influence the Cmetrics

% for PATH([1,3],:) interested in D_{3,6}
D_36 = D(:,(3-1)*c+1:3*c,6);
% for PATH([2,4],:) interested in D_{4,6}
D_46 = D(:,(4-1)*c+1:4*c,6);

pmD6 = sym(zeros(c*n,c*n)); 
pbA6D3 = sym(zeros(c*n,c*n,height(PATH))); % path specified
pA6D3 = pmD6; % no path specified

% formula for the contribution metrics following perturbation to D_{3,6}
pD36CP = sym(zeros(c*n,7,height(PATH)));
pD36CPtilde = pD36CP;
pD36Chabitat = pD36CP;

% sensitivity formula for perturbations to D_{3,6} life rates
D36sensCP = zeros(c*n,7,height(PATH)); % 7 is number of non zero life rates
D36sensCPtilde = D36sensCP;
D36sensChabitat = D36sensCP;
% elasticity formula for perturbations to D_{3,6} life rates
D36elasCP = sym(zeros(c*n,7,height(PATH))); % 7 is number of non zero life rates
D36elasCPtilde = D36elasCP;
D36elasChabitat = D36elasCP;
D36elasCP1 = sym(zeros(c*n,7,height(PATH))); % 7 is number of non zero life rates
D36elasCPtilde1 = D36elasCP1;
D36elasChabitat1 = D36elasCP1;
% perturbing life rates in D_36 one at a time
% constants used in all perturbations to D_36
E33 = zeros(n); E33(3,3)=1; % matrix E_{n,jj} j=3
vecE33 = reshape(E33,n*n,1);
for ll = 1:height(PATH)
    % constants that depend on which path is specified
    EPP = zeros(n); EPP(PATH(ll,7),PATH(ll,6))=1; 
    Z36 = reshape(kron(EPP,ones(c)),n*c*n*c,1); % constant Z_{3,k}for path ll
    
    % used in calculation of dCPtilde/dp & depends which path is specified
    bP6 = kron(eye(n),mP((PATH(ll,7)-1)*c+1:PATH(ll,7)*c,(PATH(ll,6)-1)*c+1:PATH(ll,6)*c,6));
    
    % compute sensitivities for all possible perturbations to D_36 
    for pp = 1:7 %number of non-zero entries in D_36
        syms d real     %act as delta 
        if pp==1 % perturbation is to D_36(1,3)
            % location of perturbation
            ploc = [0 0 d 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
            % 1% of life rate for elasticity substitutions
            deltasub = 268*0.01;
        elseif pp==2 % perturbation is to D_36(1,4)
            % location of perturbation
            ploc = [0 0 0 d 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
            % 1% of life rate for elasticity substitutions
            deltasub = 268*0.01;
        elseif pp==3    % perturbation is to D_36(1,5)
            % location of perturbation
            ploc = [0 0 0 0 d; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
            % 1% of life rate for elasticity substitutions
            deltasub = 89*0.01;
        elseif pp==4    % perturbation is to D_36(4,1)
            % location of perturbation
            ploc = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; d 0 0 0 0; 0 0 0 0 0];
            % 1% of life rate for elasticity substitutions
            deltasub = 0.0416*0.01;
        elseif pp==5    % perturbation is to D_36(5,3)
            % location of perturbation
            ploc = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 d 0 0];
            % 1% of life rate for elasticity substitutions
            deltasub = 0.56*0.01;
        elseif pp==6    % perturbation is to D_36(5,4)
            % location of perturbation
            ploc = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 d 0];
            % 1% of life rate for elasticity substitutions
            deltasub = 0.56*0.01;
        elseif pp==7    % perturbation is to D_36(5,5)
            % location of perturbation
            ploc = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 d];
            % 1% of life rate for elasticity substitutions
            deltasub = 0.17*0.01;
        end
        
        pD_36 = D_36 + ploc; %perturbed version of D_36
        vecDp_36 = reshape(pD_36,c*c,1);
        
        % perturbed mD_6
        for jj = 1:n 
            if jj ==3 
                pmD6((jj-1)*c+1:jj*c,(jj-1)*c+1:jj*c) = pD_36;
            else
                pmD6((jj-1)*c+1:jj*c,(jj-1)*c+1:jj*c) = mD((jj-1)*c+1:jj*c,(jj-1)*c+1:jj*c,6);
            end
        end
        
        dvecA6D3 = (Z1k(:,:,6)*X1*kron(vecE33,diff(vecDp_36))).*Z36;
        dvecA6D3habitat = kron((mS(:,:,6).*mP(:,:,6))',eye(c*n))*X1*kron(vecE33,diff(vecDp_36));
            
        D36sensCP(:,pp,ll) = X4*kronYY*dvecA6D3;
        D36sensCPtilde(:,pp,ll) = X4*kronUU*kron(bP6',eye(c*n))*dvecA6D3;
        D36sensChabitat(:,pp,ll) = X4*kronYY*dvecA6D3habitat;
        
        % perturbed contribution metrics
        pbA6D3(:,:,ll) = (pmD6(:,:)*mS(:,:,6)).*kron(EPP,ones(c)); % path specified
        pA6D3(:,:) = pmD6(:,:)*mM(:,:,6); % no path specified
        
        pD36CP(:,pp,ll) = (ones(1,c*n)*Y127*pbA6D3(:,:,ll)*Y51)';
        pD36CPtilde(:,pp,ll) = (ones(1,c*n)*Y127*pbA6D3(:,:,ll)*bP6*Y51)';
        pD36Chabitat(:,pp,ll) = (ones(1,c*n)*Y127*pA6D3(:,:)*Y51)';
        
        for ii = 1:c*n 
            % all contribution metrics will be nonzero in same place 
                % corresponding to the states populated at the start of the
                % annual cycle
            if pD36CP(ii,pp,ll) == 0 
                D36elasCP(ii,pp,ll) = 0; % or NaN
            else
                D36elasCP(ii,pp,ll) = (d/pD36CP(ii,pp,ll))*D36sensCP(ii,pp,ll);
            end
            
            if pD36CPtilde(ii,pp,ll) == 0
                D36elasCPtilde(ii,pp,ll) = 0;
            else
                D36elasCPtilde(ii,pp,ll) = (d/pD36CPtilde(ii,pp,ll))*D36sensCPtilde(ii,pp,ll);
            end
            
            if pD36Chabitat(ii,pp,ll) == 0
                D36elasChabitat(ii,pp,ll) = 0;
            else
                D36elasChabitat(ii,pp,ll) = (d/pD36Chabitat(ii,pp,ll))*D36sensChabitat(ii,pp,ll);
            end
        end
        
%         % elasticities from vectors 
%         D36elasCP1(:,pp,ll) = diag(1/pD36CP(:,pp,ll))*D36sensCP(:,pp,ll)*d;
%         D36elasCPtilde1(:,pp,ll) = diag(1/pD36CPtilde(:,pp,ll))*D36sensCPtilde(:,pp,ll)*d;
%         D36elasChabitat1(:,pp,ll) = diag(1/pD36Chabitat(:,pp,ll))*D36sensChabitat(:,pp,ll)*d;

%          % Substituting all \delta for 0.01 in elasticities function
%          D36SUBSelasCP = subs(D36elasCP,d,0.01);
%          D36SUBSelasCPtilde = subs(D36elasCPtilde,d,0.01);
%          D36SUBSelasChabitat = subs(D36elasChabitat,d,0.01); 
         
         % Substituting each \delta to be 1% of the life rate it is
         % perturbing
         D36SUBSelasCP2 = subs(D36elasCP,d,deltasub);
         D36SUBSelasCPtilde2 = subs(D36elasCPtilde,d,deltasub);
         D36SUBSelasChabitat2 = subs(D36elasChabitat,d,deltasub); 
    end            
end

%%
pmD6 = sym(zeros(c*n,c*n,height(PATH))); 
pbA6D4 = pmD6; % path specified
pA6D4 = pmD6; % no path specified
% formula for the contribution metrics following perturbation to D_{4,6}
pD46CP = sym(zeros(c*n,height(PATH)));
pD46CPtilde = pD46CP;
pD46Chabitat = pD46CP;
% sensitivity formula for perturbations to D_{4,6} life rates
D46sensCP = sym(zeros(c*n,7,height(PATH)));
D46sensCPtilde = D46sensCP;
D46sensChabitat = D46sensCP; 
% elasticity formula for perturbations to D_{4,6} life rates
D46elasCP = sym(zeros(c*n,7,height(PATH))); % 7 is number of non zero life rates
D46elasCPtilde = D46elasCP;
D46elasChabitat = D46elasCP;
% perturbing life rates in D_46 one at a time
% constants used in all perturbations to D_46
E44 = zeros(n); E44(4,4)=1; % matrix E_{n,jj} j=4
vecE44 = reshape(E44,n*n,1);
for ll = 1:height(PATH)
    % constants that depend on which path is specified
    EPP = zeros(n); EPP(PATH(ll,7),PATH(ll,6))=1; 
    Z36 = reshape(kron(EPP,ones(c)),n*c*n*c,1); % constant Z_{3,k}for path ll
    
    % used in calculation of dCPtilde/dp & depends which path is specified
    bP6 = kron(eye(n),mP((PATH(ll,7)-1)*c+1:PATH(ll,7)*c,(PATH(ll,6)-1)*c+1:PATH(ll,6)*c,6));
    
    % compute sensitivities for all possible perturbations to D_46 
    for pp = 1:7 %number of non-zero entries in D_46
        syms d real     %act as delta 
        if pp==1 % perturbation is to D_36(1,3)
            % location of perturbation
            ploc = [0 0 d 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
            % 1% of life rate for elasticity substitutions
            deltasub = 268*0.01;
        elseif pp==2 % perturbation is to D_36(1,4)
            % location of perturbation
            ploc = [0 0 0 d 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
            % 1% of life rate for elasticity substitutions
            deltasub = 268*0.01;
        elseif pp==3    % perturbation is to D_36(1,5)
            % location of perturbation
            ploc = [0 0 0 0 d; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
            % 1% of life rate for elasticity substitutions
            deltasub = 89*0.01;
        elseif pp==4    % perturbation is to D_36(4,1)
            % location of perturbation
            ploc = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; d 0 0 0 0; 0 0 0 0 0];
            % 1% of life rate for elasticity substitutions
            deltasub = 0.0416*0.01;
        elseif pp==5    % perturbation is to D_36(5,3)
            % location of perturbation
            ploc = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 d 0 0];
            % 1% of life rate for elasticity substitutions
            deltasub = 0.56*0.01;
        elseif pp==6    % perturbation is to D_36(5,4)
            % location of perturbation
            ploc = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 d 0];
            % 1% of life rate for elasticity substitutions
            deltasub = 0.56*0.01;
        elseif pp==7    % perturbation is to D_36(5,5)
            % location of perturbation
            ploc = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 d];
            % 1% of life rate for elasticity substitutions
            deltasub = 0.17*0.01;
        end
        
        pD_46 = D_46 + ploc; %perturbed version of D_46
        vecDp_46 = reshape(pD_46,c*c,1);
        
        % perturbed mD_6
        for jj = 1:n 
            if jj ==4 
                pmD6((jj-1)*c+1:jj*c,(jj-1)*c+1:jj*c,ll) = pD_46;
            else
                pmD6((jj-1)*c+1:jj*c,(jj-1)*c+1:jj*c,ll) = mD((jj-1)*c+1:jj*c,(jj-1)*c+1:jj*c,6);
            end
        end
        
        dvecA6D4 = (Z1k(:,:,6)*X1*kron(vecE44,diff(vecDp_46))).*Z36;
        dvecA6D4habitat = kron((mS(:,:,6).*mP(:,:,6))',eye(c*n))*X1*kron(vecE44,diff(vecDp_46));
            
        D46sensCP(:,pp,ll) = X4*kronYY*dvecA6D4;
        D46sensCPtilde(:,pp,ll) = X4*kronUU*kron(bP6',eye(c*n))*dvecA6D4;
        D46sensChabitat(:,pp,ll) = X4*kronYY*dvecA6D4habitat;
        
        % perturbed contribution metrics
        pbA6D4(:,:,ll) = (pmD6(:,:,ll)*mS(:,:,6)).*kron(EPP,ones(c)); % path specified
        pA6D4(:,:,ll) = pmD6(:,:,ll)*mM(:,:,6); % no path specified
        
        pD46CP(:,ll) = (ones(1,c*n)*Y127*pbA6D4(:,:,ll)*Y51)';
        pD46CPtilde(:,ll) = (ones(1,c*n)*Y127*pbA6D4(:,:,ll)*bP6*Y51)';
        pD46Chabitat(:,ll) = (ones(1,c*n)*Y127*pA6D4(:,:,ll)*Y51)';
        
        for ii = 1:c*n 
            % all contribution metrics will be nonzero in same place 
                % corresponding to the states populated at the start of the
                % annual cycle
            if pD46CP(ii,ll) == 0 
                D46elasCP(ii,pp,ll) = 0; % or NaN
            else
                D46elasCP(ii,pp,ll) = (d/pD46CP(ii,ll))*D46sensCP(ii,pp,ll);
            end
            
            if pD46CPtilde(ii,ll) == 0
                D46elasCPtilde(ii,pp,ll) = 0;
            else
                D46elasCPtilde(ii,pp,ll) = (d/pD46CPtilde(ii,ll))*D46sensCPtilde(ii,pp,ll);
            end
            
            if pD46Chabitat(ii,ll) == 0
                D46elasChabitat(ii,pp,ll) = 0;
            else
                D46elasChabitat(ii,pp,ll) = (d/pD46Chabitat(ii,ll))*D46sensChabitat(ii,pp,ll);
            end
        end
        
%         % Substituting all \delta for 0.01 in elasticities function
%          D46SUBSelasCP = subs(D46elasCP,d,0.01);
%          D46SUBSelasCPtilde = subs(D46elasCPtilde,d,0.01);
%          D46SUBSelasChabitat = subs(D46elasChabitat,d,0.01);  

        % Substituting each \delta to be 1% of the life rate it is
         % perturbing
         D46SUBSelasCP2 = subs(D46elasCP,d,deltasub);
         D46SUBSelasCPtilde2 = subs(D46elasCPtilde,d,deltasub);
         D46SUBSelasChabitat2 = subs(D46elasChabitat,d,deltasub); 
    end            
end

%% SENSITIVITY OF CONTRIBUTION METRICS WITH RESPECT TO MOVEMENT SURVIVAL RATES
% MIGRATION HAPPENS BEFORE DEMOGRAPHY
% assume stages 1, 2 and 3 cannnot alter their movement survival rates
% there are 3 life rates that could change in both stage 4 and 5, so we
% will do a total of 6 perturbations

pmS46 = sym(zeros(c*n,c*n,height(PATH))); 
pmS56 = pmS46; 
%
pbA6S4 = pmS46; % path specified
pA6S4 = pmS46; % path not specified
pbA6S5 = pmS56; % path specified
pA6S5 = pmS56; % path not specified
% formula for the contribution metrics following perturbation to S^4_6
pS46CP = sym(zeros(c*n,height(PATH)));
pS46CPtilde = pS46CP;
pS46Chabitat = pS46CP;
% formula for the contribution metrics following perturbation to S^5_6
pS56CP = sym(zeros(c*n,height(PATH)));
pS56CPtilde = pS56CP;
pS56Chabitat = pS56CP;

% sensitivity formula for perturbations to S^4_6 
S46sensCP = zeros(c*n,3,height(PATH)); 
S46sensCPtilde = S46sensCP;
S46sensChabitat = S46sensCP;
% sensitivity formula for perturbations to S^5_6 
S56sensCP = zeros(c*n,3,height(PATH));
S56sensCPtilde = S56sensCP;
S56sensChabitat = S56sensCP; 
% elasticity formula for perturbations to S^4_6 life rates
S46elasCP = sym(zeros(c*n,3,height(PATH))); % 3 is number of non-zero life rates
S46elasCPtilde = S46elasCP;
S46elasChabitat = S46elasCP;
% elasticity formula for perturbations to S^5_6 life rates
S56elasCP = sym(zeros(c*n,3,height(PATH))); 
S56elasCPtilde = S56elasCP;
S56elasChabitat = S56elasCP;
% perturbing movement rates associated with specified paths if non-zero
%%% assumption that new paths cannot be created 
for ll = 1:height(PATH)
    % constants that depend on which path is specified
    EPP = zeros(n); EPP(PATH(ll,7),PATH(ll,6))=1; 
    Z36 = reshape(kron(EPP,ones(c)),n*c*n*c,1);%constant Z_{3,k}for path ll
    
    % used in calculation of dCPtilde/dp & depends which path is specified
    bP6 = kron(eye(n),mP((PATH(ll,7)-1)*c+1:PATH(ll,7)*c,(PATH(ll,6)-1)*c+1:PATH(ll,6)*c,6));
    
    for pp = 1:6
        syms d real % act as delta
        if pp <= 3 % first 3 perturbations will be to stage 4 
            Ecii = zeros(c); Ecii(4,4)=1; %matrix E_{c,ii} i=4
            vecEcii = reshape(Ecii,c*c,1); 
            
            if pp == 1  % perturbation is to S^4_6(3,2)
                ploc = [0 0 0 0; 0 0 0 0; 0 d 0 0; 0 0 0 0];
                % 1% of life rate for elasticity substitutions
                deltasub = 0.733*0.01;
            elseif pp == 2 % perturbation is to S^4_6(4,2)
                ploc = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 d 0 0];
                % 1% of life rate for elasticity substitutions
                deltasub = 0.544*0.01;
            elseif pp == 3 % perturbation is to S^4_6(4,3)
                ploc = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 d 0];
                % 1% of life rate for elasticity substitutions
                deltasub = 0.742*0.01;
            end 
            
            pS46 = S(:,3*n+1:4*n,6) + ploc; % perturbed version of S^4_6
            vecSp46 = reshape(pS46,n*n,1); 
            
            % perturbed mS_6
            for ii = 1:c 
                E = zeros(c); E(ii,ii) = 1;  
                if ii == 4 
                    pmS46(:,:,ll) = pmS46(:,:,ll) + kron(pS46,E);
                else
                    pmS46(:,:,ll) = pmS46(:,:,ll) + kron(S(:,(ii-1)*n+1:ii*n,6),E);
                end
            end
            
            dvecA6S4 = (X2k(:,:,6)*((X1*kron(diff(vecSp46),vecEcii)).*vecmP6)).*Z36;
            dvecA6S4habitat = X2k(:,:,6)*(X1*kron(diff(vecSp46),vecEcii).*vecmP6);
            
            S46sensCP(:,pp,ll) = X4*kronYY*dvecA6S4;
            S46sensCPtilde(:,pp,ll) = X4*kronUU*kron(bP6',eye(c*n))*dvecA6S4;
            S46sensChabitat(:,pp,ll) = X4*kronYY*dvecA6S4habitat;
            
            % perturbed contribution metrics
            pbA6S4(:,:,ll) = (mD(:,:,6)*pmS46(:,:,ll)).*kron(EPP,ones(c)); % path specified
            pA6D4(:,:,ll) = mD(:,:,6)*(pmS46(:,:,ll).*mP(:,:,6)); % no path specified

            pS46CP(:,ll) = (ones(1,c*n)*Y127*pbA6S4(:,:,ll)*Y51)';
            pS46CPtilde(:,ll) = (ones(1,c*n)*Y127*pbA6S4(:,:,ll)*bP6*Y51)';
            pS46Chabitat(:,ll) = (ones(1,c*n)*Y127*pA6S4(:,:,ll)*Y51)';
            
            for ii = 1:c*n 
                % all contribution metrics will be nonzero in same place 
                    % corresponding to the states populated at the start of the
                    % annual cycle
                if pS46CP(ii,ll) == 0 
                    S46elasCP(ii,pp,ll) = 0; % or NaN
                else
                    S46elasCP(ii,pp,ll) = (d/pS46CP(ii,ll))*S46sensCP(ii,pp,ll);
                end

                if pS46CPtilde(ii,ll) == 0
                    S46elasCPtilde(ii,pp,ll) = 0;
                else
                    S46elasCPtilde(ii,pp,ll) = (d/pS46CPtilde(ii,ll))*S46sensCPtilde(ii,pp,ll);
                end

                if pS46Chabitat(ii,ll) == 0
                    S46elasChabitat(ii,pp,ll) = 0;
                else
                    S46elasChabitat(ii,pp,ll) = (d/pS46Chabitat(ii,ll))*S46sensChabitat(ii,pp,ll);
                end
            end
              
        else % last 3 perturbations will be to stage 5
            Ecii = zeros(c); Ecii(5,5)=1; %matrix E_{c,ii} i=5
            vecEcii = reshape(Ecii,c*c,1); 
            
            if pp == 4  % perturbation is to S^5_6(3,2)
                ploc = [0 0 0 0; 0 0 0 0; 0 d 0 0; 0 0 0 0];
                % 1% of life rate for elasticity substitutions
                deltasub = 0.733*0.01;
            elseif pp == 5 % perturbation is to S^5_6(4,2)
                ploc = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 d 0 0];
                % 1% of life rate for elasticity substitutions
                deltasub = 0.544*0.01;
            elseif pp == 6 % perturbation is to S^5_6(4,3)
                ploc = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 d 0];
                % 1% of life rate for elasticity substitutions
                deltasub = 0.742*0.01;
            end 
            
            pS56 = S(:,4*n+1:5*n,6) + ploc; % perturbed version of S^5_6
            vecSp56 = reshape(pS56,n*n,1); 
            
            % perturbed mS_6
            for ii = 1:c 
                E = zeros(c); E(ii,ii) = 1;  
                if ii == 5 
                    pmS56(:,:,ll) = pmS56(:,:,ll) + kron(pS56,E);
                else
                    pmS56(:,:,ll) = pmS56(:,:,ll) + kron(S(:,(ii-1)*n+1:ii*n,6),E);
                end
            end
            
            dvecA6S5 = (X2k(:,:,6)*((X1*kron(diff(vecSp56),vecEcii)).*vecmP6)).*Z36;
            dvecA6S5habitat = X2k(:,:,6)*(X1*kron(diff(vecSp56),vecEcii).*vecmP6);
            
            S56sensCP(:,pp-3,ll) = X4*kronYY*dvecA6S5;
            S56sensCPtilde(:,pp-3,ll) = X4*kronUU*kron(bP6',eye(c*n))*dvecA6S5;
            S56sensChabitat(:,pp-3,ll) = X4*kronYY*dvecA6S5habitat;
            
            % perturbed contribution metrics
            pbA6S5(:,:,ll) = (mD(:,:,6)*pmS56(:,:,ll)).*kron(EPP,ones(c)); % path specified
            pA6S5(:,:,ll) = mD(:,:,6)*(pmS56(:,:,ll).*mP(:,:,6)); % no path specified

            pS56CP(:,ll) = (ones(1,c*n)*Y127*pbA6S5(:,:,ll)*Y51)';
            pS56CPtilde(:,ll) = (ones(1,c*n)*Y127*pbA6S5(:,:,ll)*bP6*Y51)';
            pS56Chabitat(:,ll) = (ones(1,c*n)*Y127*pA6S5(:,:,ll)*Y51)';
            
            for ii = 1:c*n 
                % all contribution metrics will be nonzero in same place 
                    % corresponding to the states populated at the start of the
                    % annual cycle
                if pS56CP(ii,ll) == 0 
                    S56elasCP(ii,pp-3,ll) = 0; % or NaN
                else
                    S56elasCP(ii,pp-3,ll) = (d/pS56CP(ii,ll))*S56sensCP(ii,pp-3,ll);
                end

                if pS56CPtilde(ii,ll) == 0
                    S56elasCPtilde(ii,pp-3,ll) = 0;
                else
                    S56elasCPtilde(ii,pp-3,ll) = (d/pS56CPtilde(ii,ll))*S56sensCPtilde(ii,pp-3,ll);
                end

                if pS56Chabitat(ii,ll) == 0
                    S56elasChabitat(ii,pp-3,ll) = 0;
                else
                    S56elasChabitat(ii,pp-3,ll) = (d/pS56Chabitat(ii,ll))*S56sensChabitat(ii,pp-3,ll);
                end
            end
        end
        
%         % Substituting all \delta for 0.01 in elasticities function
%          S46SUBSelasCP = subs(S46elasCP,d,0.01);
%          S46SUBSelasCPtilde = subs(S46elasCPtilde,d,0.01);
%          S46SUBSelasChabitat = subs(S46elasChabitat,d,0.01);  
%          %
%          S56SUBSelasCP = subs(S56elasCP,d,0.01);
%          S56SUBSelasCPtilde = subs(S56elasCPtilde,d,0.01);
%          S56SUBSelasChabitat = subs(S56elasChabitat,d,0.01);  
         
         % Substituting each \delta to be 1% of the life rate it is
         % perturbing
         S46SUBSelasCP2 = subs(S46elasCP,d,deltasub);
         S46SUBSelasCPtilde2 = subs(S46elasCPtilde,d,deltasub);
         S46SUBSelasChabitat2 = subs(S46elasChabitat,d,deltasub); 
         %
         S56SUBSelasCP2 = subs(S56elasCP,d,deltasub);
         S56SUBSelasCPtilde2 = subs(S56elasCPtilde,d,deltasub);
         S56SUBSelasChabitat2 = subs(S56elasChabitat,d,deltasub); 
    end 
end


%% SENSITIVITY OF CONTRIBUTION METRICS WITH RESPECT TO MOVEMENT SURVIVAL RATES V2

%%%% UPDATE ELASTICITY EQUATIONS REQUIRED !!! 


% MIGRATION HAPPENS BEFORE DEMOGRAPHY
% assume stages 1, 2 and 3 cannnot alter their movement survival rates
% there are 3 life rates that could change in both stage 4 and 5.
% HERE we assume that the perturbation to the movement survival affects
% stage 4 and 5 in the same way. Motivated by the unperturbed movement
% survival rates being the same for stages 4 and 5. 

pmS = sym(zeros(c*n,c*n,height(PATH))); 
%
pbA6S = pmS; % path specified
pA6S = pmS; % path not specified
% formula for the contribution metrics following perturbation to S^4_6 and S^5_6
pSCP = sym(zeros(c*n,height(PATH)));
pSCPtilde = pSCP;
pSChabitat = pSCP;

% sensitivity formula for perturbations to S^4_6 and S^5_6
SsensCP = zeros(c*n,3,height(PATH)); 
SsensCPtilde = SsensCP;
SsensChabitat = SsensCP;
% elasticity formula for perturbations to S^4_6 and S^5_6 life rates
SelasCP = sym(zeros(c*n,3,height(PATH))); 
SelasCPtilde = SelasCP;
SelasChabitat = SelasCP;
% perturbing movement rates associated with specified paths if non-zero
%%% assumption that new paths cannot be created 
for ll = 1:height(PATH)
    % constants that depend on which path is specified
    EPP = zeros(n); EPP(PATH(ll,7),PATH(ll,6))=1; 
    Z36 = reshape(kron(EPP,ones(c)),n*c*n*c,1);%constant Z_{3,k}for path ll
    
    % used in calculation of dCPtilde/dp & depends which path is specified
    bP6 = kron(eye(n),mP((PATH(ll,7)-1)*c+1:PATH(ll,7)*c,(PATH(ll,6)-1)*c+1:PATH(ll,6)*c,6));
    
    for pp = 1:3
        syms d real% act as delta
              
        if pp == 1  % perturbation is to S^4_6(3,2) and S^5_6(3,2)
            ploc = [0 0 0 0; 0 0 0 0; 0 d 0 0; 0 0 0 0];
        elseif pp == 2 % perturbation is to S^4_6(4,2) and S^5_6(4,2)
            ploc = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 d 0 0];
        elseif pp == 3 % perturbation is to S^4_6(4,3) and S^5_6(4,3)
            ploc = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 d 0];
        end 

        pS46 = S(:,3*n+1:4*n,6) + ploc; % perturbed version of S^4_6
        vecSp46 = reshape(pS46,n*n,1); 
        
        pS56 = S(:,4*n+1:5*n,6) + ploc; % perturbed version of S^5_6
        vecSp56 = reshape(pS56,n*n,1);
        
        % perturbed mS_6
        for ii = 1:c 
            E = zeros(c); E(ii,ii) = 1;  
            if ii == 4  % stage 4 perturbed
                pmS(:,:,ll) = pmS(:,:,ll) + kron(pS46,E);
            elseif ii == 5 % stage 5 perturbed
                pmS(:,:,ll) = pmS(:,:,ll) + kron(pS56,E);
            else
                pmS(:,:,ll) = pmS(:,:,ll) + kron(S(:,(ii-1)*n+1:ii*n,6),E);
            end
        end
        
        % parts within the sum over i. X1 is constant so can be pulled out
        kronpS4 = kron(diff(vecSp46),vecEc44);
        kronpS5 = kron(diff(vecSp56),vecEc55);
        
        dvecA6S = (X2k(:,:,6)*((X1*(kronpS4+kronpS5)).*vecmP6)).*Z36;
        dvecA6Shabitat = X2k(:,:,6)*(X1*(kronpS4+kronpS5).*vecmP6);

        SsensCP(:,pp,ll) = X4*kronYY*dvecA6S;
        SsensCPtilde(:,pp,ll) = X4*kronUU*kron(bP6',eye(c*n))*dvecA6S;
        SsensChabitat(:,pp,ll) = X4*kronYY*dvecA6Shabitat;
        
        % perturbed contribution metrics
        pbA6S(:,:,ll) = (mD(:,:,6)*pmS(:,:,ll)).*kron(EPP,ones(c)); % path specified
        pA6S(:,:,ll) = mD(:,:,6)*(pmS(:,:,ll).*mP(:,:,6)); % no path specified

        pSCP(:,ll) = (ones(1,c*n)*Y127*pbA6S(:,:,ll)*Y51)';
        pSCPtilde(:,ll) = (ones(1,c*n)*Y127*pbA6S(:,:,ll)*bP6*Y51)';
        pSChabitat(:,ll) = (ones(1,c*n)*Y127*pA6S(:,:,ll)*Y51)';

        for ii = 1:c*n 
            % all contribution metrics will be nonzero in same place 
                % corresponding to the states populated at the start of the
                % annual cycle
            if pSCP(ii,ll) == 0 
                SelasCP(ii,pp,ll) = 0; % or NaN
            else
                SelasCP(ii,pp,ll) = (d/pSCP(ii,ll))*SsensCP(ii,pp,ll);
            end

            if pSCPtilde(ii,ll) == 0
                SelasCPtilde(ii,pp,ll) = 0;
            else
                SelasCPtilde(ii,pp,ll) = (d/pSCPtilde(ii,ll))*SsensCPtilde(ii,pp,ll);
            end

            if pSChabitat(ii,ll) == 0
                SelasChabitat(ii,pp,ll) = 0;
            else
                SelasChabitat(ii,pp,ll) = (d/pSChabitat(ii,ll))*SsensChabitat(ii,pp,ll);
            end
        end
        
        % Substituting all \delta for 0.01 in elasticities function
        SSUBSelasCP = subs(SelasCP,d,0.01);
        SSUBSelasCPtilde = subs(SelasCPtilde,d,0.01);
        SSUBSelasChabitat = subs(SelasChabitat,d,0.01);  
    end
end

% this gives the same result as summing S46sensCP and S56sensCP, due to X1
% being constant and so can pull outside of the sum.

%% SENSITIVITY OF CONTRIBUTION METRICS WITH RESPECT TO PROPORTION MIGRATING RATES
% MIGRATION HAPPENS BEFORE DEMOGRAPHY 
% there are 4 proportion rates for both stage 4 & 5, but within each stage,
% the proportion rates can only be perturbed in pairs, as the columns of
% P^i_k need to sum to 1. Thus there are only two available perturbations
% for each stage, and so 4 perturbations in total. 

% the sensitivity of \bC(\sP) is 0 for all as there is a specified path in
% season 6 (Z_{2,6}=0). season 6 is also when we are performing the
% perturbations. 


pmP46 = sym(zeros(c*n,c*n,height(PATH))); 
pmP56 = pmP46; 
% bold P terms are affected by perturbation to P^i_k matrices
pbP6P4 = pmP46; % path specified 
pbP6P5 = pbP6P4; % path specified
% formula for the contribution metrics following perturbation to P^4_6
pP46CP = sym(zeros(c*n,height(PATH)));
pP46CPtilde = pP46CP;
pP46Chabitat = pP46CP;
% formula for the contribution metrics following perturbation to P^5_6
pP56CP = sym(zeros(c*n,height(PATH)));
pP56CPtilde = pP56CP;
pP56Chabitat = pP56CP;

% sensitivity formula for perturbations to P^4_6
P46sensCPtilde = sym(zeros(c*n,2,height(PATH))); 
P46sensChabitat = P46sensCPtilde;
% elasticity formula for perturbations to P^4_6 life rates
P46elasCPtilde = P46sensCPtilde;
P46elasChabitat = P46elasCPtilde;
% sensitivity formula for perturbations to P^5_6 
P56sensCPtilde = sym(zeros(c*n,2,height(PATH))); 
P56sensChabitat = P56sensCPtilde;
% elasticity formula for perturbations to P^5_6 life rates
P56elasCPtilde = P56sensCPtilde;
P56elasChabitat = P56elasCPtilde;
% dvec(A_k)/dp^T is zero as does not contain P matrices 

for ll = 1:height(PATH)
    % constants that depend on which path is specified   
    EPP = zeros(n); EPP(PATH(ll,7),PATH(ll,6))=1; 
    bA6 = Ahat(:,:,6).*kron(EPP,ones(c));
    
    T_16 = zeros(c,c*n); T_16(:,(PATH(ll,7)-1)*c+1:PATH(ll,7)*c) = eye(c);
    T_26 = zeros(c*n,c); T_26((PATH(ll,6)-1)*c+1:PATH(ll,6)*c,:) = eye(c);
    
    for pp = 1:2
        syms d real % act as delta
        
        if pp == 1 
    % perturbation is to P^4_6(3,2) and P^4_6(4,2) and stage 5 equivalents
            ploc = [0 0 0 0; 0 0 0 0; 0 d 0 0; 0 -d 0 0];
            % 1% of life rate being positively perturbed for elasticity substitutions
            deltasub = 0.556*0.01;
        else
    % perturbation is to P^4_6(3,3) and P^4_6(4,3) and stage 5 equivalents
            ploc = [0 0 0 0; 0 0 0 0; 0 0 d 0; 0 0 -d 0];
            % 1% of life rate being positively perturbed for elasticity substitutions
            deltasub = 0.532*0.01;
        end
        
        pP46 = P(:,3*n+1:4*n,6) + ploc; % perturbed version of P^4_6
        vecPp46 = reshape(pP46,n*n,1);
        
        pP56 = P(:,4*n+1:5*n,6) + ploc; % perturbed version of P^5_6
        vecPp56 = reshape(pP56,n*n,1);
        
        % perturbed mP_6
        for ii = 1:c 
            E = zeros(c); E(ii,ii) = 1;  
            if ii == 4 
                pmP46(:,:,ll) = pmP46(:,:,ll) + kron(pP46,E);
                pmP56(:,:,ll) = pmP56(:,:,ll) + kron(P(:,(ii-1)*n+1:ii*n,6),E);
            elseif ii == 5 
                pmP46(:,:,ll) = pmP46(:,:,ll) + kron(P(:,(ii-1)*n+1:ii*n,6),E);
                pmP56(:,:,ll) = pmP56(:,:,ll) + kron(pP56,E);
            else
                pmP46(:,:,ll) = pmP46(:,:,ll) + kron(P(:,(ii-1)*n+1:ii*n,6),E);
                pmP56(:,:,ll) = pmP56(:,:,ll) + kron(P(:,(ii-1)*n+1:ii*n,6),E);
            end
        end
        
        P46sensCPtilde(:,pp,ll) = X4*kronUU*kron(eye(c*n),bA6)*(Z46*X1*...
            kron(vecIn,(kron(T_26',T_16)*(X1*kron(diff(vecPp46),vecEc44)))));
        P46sensChabitat(:,pp,ll) = X4*kronYY*X2k(:,:,6)*(vecmS6.*(X1*kron(diff(vecPp46),vecEc44)));
        
        P56sensCPtilde(:,pp,ll) = X4*kronUU*kron(eye(c*n),bA6)*(Z46*X1*...
            kron(vecIn,(kron(T_26',T_16)*(X1*kron(diff(vecPp56),vecEc55)))));
        P56sensChabitat(:,pp,ll) = X4*kronYY*X2k(:,:,6)*(vecmS6.*(X1*kron(diff(vecPp56),vecEc55)));
        
        % perturbed contribution metrics
        pbP6P4(:,:,ll) = kron(eye(n),T_16*pmP46(:,:,ll)*T_26);
        pbP6P5(:,:,ll) = kron(eye(n),T_16*pmP56(:,:,ll)*T_26);
        
        pA6P4 = mD(:,:,6)*(mS(:,:,6).*pmP46(:,:,ll));
        pA6P5 = mD(:,:,6)*(mS(:,:,6).*pmP56(:,:,ll));
        
        pP46CP(:,ll) = (ones(1,c*n)*Y127*bA6*Y51)';
        pP46CPtilde(:,ll) = (ones(1,c*n)*Y127*bA6*pbP6P4(:,:,ll)*Y51)';
        pP46Chabitat(:,ll) = (ones(1,c*n)*Y127*pA6P4*Y51)';

        pP56CP(:,ll) = (ones(1,c*n)*Y127*bA6*Y51)';
        pP56CPtilde(:,ll) = (ones(1,c*n)*Y127*bA6*pbP6P5(:,:,ll)*Y51)';
        pP56Chabitat(:,ll) = (ones(1,c*n)*Y127*pA6P5*Y51)';
        
        for ii = 1:c*n 
                % all contribution metrics will be nonzero in same place 
                    % corresponding to the states populated at the start of the
                    % annual cycle
                if pP46CPtilde(ii,ll) == 0
                    P46elasCPtilde(ii,pp,ll) = 0;
                else
                    P46elasCPtilde(ii,pp,ll) = (d/pP46CPtilde(ii,ll))*P46sensCPtilde(ii,pp,ll);
                end

                if pP46Chabitat(ii,ll) == 0
                    P46elasChabitat(ii,pp,ll) = 0;
                else
                    P46elasChabitat(ii,pp,ll) = (d/pP46Chabitat(ii,ll))*P46sensChabitat(ii,pp,ll);
                end
                
                if pP56CPtilde(ii,ll) == 0
                    P56elasCPtilde(ii,pp,ll) = 0;
                else
                    P56elasCPtilde(ii,pp,ll) = (d/pP56CPtilde(ii,ll))*P56sensCPtilde(ii,pp,ll);
                end

                if pP56Chabitat(ii,ll) == 0
                    P56elasChabitat(ii,pp,ll) = 0;
                else
                    P56elasChabitat(ii,pp,ll) = (d/pP56Chabitat(ii,ll))*P56sensChabitat(ii,pp,ll);
                end
        end
%         % Substituting all \delta for 0.01 in elasticities function
%          P46SUBSelasCP = subs(P46elasCP,d,0.01);
%          P46SUBSelasCPtilde = subs(P46elasCPtilde,d,0.01);
%          P46SUBSelasChabitat = subs(P46elasChabitat,d,0.01);  
%          %
%          P56SUBSelasCP = subs(P56elasCP,d,0.01);
%          P56SUBSelasCPtilde = subs(P56elasCPtilde,d,0.01);
%          P56SUBSelasChabitat = subs(P56elasChabitat,d,0.01);  
         
         % Substituting each \delta to be 1% of the life rate it is
         % being positively perturbed
%          P46SUBSelasCP2 = subs(P46elasCP,d,deltasub);
         P46SUBSelasCPtilde2 = subs(P46elasCPtilde,d,deltasub);
         P46SUBSelasChabitat2 = subs(P46elasChabitat,d,deltasub); 
         %
%          P56SUBSelasCP2 = subs(P56elasCP,d,deltasub);
         P56SUBSelasCPtilde2 = subs(P56elasCPtilde,d,deltasub);
         P56SUBSelasChabitat2 = subs(P56elasChabitat,d,deltasub); 
    end
end 

%% SENSITIVITY OF CONTRIBUTION METRICS WITH RESPECT TO PROPORTION MIGRATING V2

%%%% STILL NEEDS SOME FIGURING OUT. habitat contribution metrics appear to
%%%% give correct  result for sensitivities. Using T_{1,6} =
%%%% kron(ones(1,n),eye(c)) does not maintain the structure required. So
%%%% seems as if summing sensitivities of  Ctilde once calculated is the
%%%% correct way. Elasticities of habitat contribution metrics does not
%%%% give same result (or magnitude) as summing the elasticities of Ctilde,
%%%% but sign of elasticities of habitat metrics is the same as the sign of
%%%% the corresponding sensitivities -- this is not the case for summing
%%%% the elasticities of the Ctildes.

 % NOW ONLY STARTING HABITAT OF SEASON 6 IS SPECIFIED 
PATHstart = zeros(2,s+1); 
PATHstart(1,6) = 2; % individuals that start season 6 in habitat 2
PATHstart(2,6) = 3; % individuals that start season 6 in habitat 3

pmP46start = sym(zeros(c*n,c*n,height(PATHstart))); 
pmP56start = pmP46start; 
% bold P terms are affected by perturbation to P^i_k matrices
pbP6P4start = pmP46start; % path specified 
pbP6P5start = pbP6P4start; % path specified
% formula for the contribution metrics following perturbation to P^4_6
pP46CPtildestart =sym(zeros(c*n,height(PATHstart)));
pP46Chabitatstart = pP46CPtildestart;
% formula for the contribution metrics following perturbation to P^5_6
pP56CPtildestart = sym(zeros(c*n,height(PATHstart)));
pP56Chabitatstart = pP56CPtildestart;

% sensitivity formula for perturbations to P^4_6
P46sensCPtildestart = zeros(c*n,2,height(PATHstart)); 
P46sensChabitatstart = P46sensCPtildestart;
% elasticity formula for perturbations to P^4_6 life rates
P46elasCPstart = sym(zeros(c*n,2,height(PATHstart))); 
P46elasCPtildestart = P46elasCPstart;
P46elasChabitatstart = P46elasCPstart;
% sensitivity formula for perturbations to P^5_6 
P56sensCPtildestart = zeros(c*n,2,height(PATHstart)); 
P56sensChabitatstart = P56sensCPtildestart;
% elasticity formula for perturbations to P^5_6 life rates
P56elasCPstart = sym(zeros(c*n,2,height(PATHstart))); 
P56elasCPtildestart = P56elasCPstart;
P56elasChabitatstart = P56elasCPstart;
% dvec(A_k)/dp^T is zero as does not contain P matrices 

for ll = 1:height(PATHstart)
    % constants that depend on which path is specified   
    EPP = zeros(n); EPP(:,PATHstart(ll,6))=1; % can start in all habitats in season 7
    bA6 = Ahat(:,:,6).*kron(EPP,ones(c));
    
    T_16start = kron(ones(1,n),eye(c));
    e_n = zeros(n,1); e_n(PATHstart(ll,6)) = 1;
    T_26start = kron(e_n,eye(c));
    
    for pp = 1:2
        syms d real % act as delta
        
        if pp == 1 
    % perturbation is to P^4_6(3,2) and P^4_6(4,2) and stage 5 equivalents
            ploc = [0 0 0 0; 0 0 0 0; 0 d 0 0; 0 -d 0 0];
    % 1% of life rate being positively perturbed for elasticity substitutions
            deltasub = 0.556*0.01;
        else
    % perturbation is to P^4_6(3,3) and P^4_6(4,3) and stage 5 equivalents
            ploc = [0 0 0 0; 0 0 0 0; 0 0 d 0; 0 0 -d 0];
            % 1% of life rate being positively perturbed for elasticity substitutions
            deltasub = 0.532*0.01;
        end
        
        pP46start = P(:,3*n+1:4*n,6) + ploc; % perturbed version of P^4_6
        vecpP46start = reshape(pP46start,n*n,1);
        
        pP56start = P(:,4*n+1:5*n,6) + ploc; % perturbed version of P^5_6
        vecpP56start = reshape(pP56start,n*n,1);
        
        % perturbed mP_6
        for ii = 1:c 
            E = zeros(c); E(ii,ii) = 1;  
            if ii == 4 
                pmP46start(:,:,ll) = pmP46start(:,:,ll) + kron(pP46start,E);
                pmP56start(:,:,ll) = pmP56start(:,:,ll) + kron(P(:,(ii-1)*n+1:ii*n,6),E);
            elseif ii == 5 
                pmP46start(:,:,ll) = pmP46start(:,:,ll) + kron(P(:,(ii-1)*n+1:ii*n,6),E);
                pmP56start(:,:,ll) = pmP56start(:,:,ll) + kron(pP56start,E);
            else
                pmP46start(:,:,ll) = pmP46start(:,:,ll) + kron(P(:,(ii-1)*n+1:ii*n,6),E);
                pmP56start(:,:,ll) = pmP56start(:,:,ll) + kron(P(:,(ii-1)*n+1:ii*n,6),E);
            end
        end
        
        P46sensCPtildestart(:,pp,ll) = X4*kronUU*kron(eye(c*n),bA6)*(Z46*X1*...
            kron(vecIn,(kron(T_26start',T_16start)*(X1*kron(diff(vecpP46start),vecEc44)))));
        P46sensChabitatstart(:,pp,ll) = X4*kronYY*X2k(:,:,6)*(vecmS6.*(X1*kron(diff(vecpP46start),vecEc44)));
        
        P56sensCPtildestart(:,pp,ll) = X4*kronUU*kron(eye(c*n),bA6)*(Z46*X1*...
            kron(vecIn,(kron(T_26start',T_16start)*(X1*kron(diff(vecpP56start),vecEc55)))));
        P56sensChabitatstart(:,pp,ll) = X4*kronYY*X2k(:,:,6)*(vecmS6.*(X1*kron(diff(vecpP56start),vecEc55)));
        
        % perturbed contribution metrics
        pbP6P4start(:,:,ll) = kron(eye(n),T_16start*pmP46start(:,:,ll)*T_26start);
        pbP6P5start(:,:,ll) = kron(eye(n),T_16start*pmP56start(:,:,ll)*T_26start);
        
        pA6P4start = mD(:,:,6)*(mS(:,:,6).*pmP46start(:,:,ll));
        pA6P5start = mD(:,:,6)*(mS(:,:,6).*pmP56start(:,:,ll));
        
%         pP46CP(:,ll) = (ones(1,c*n)*Y127*bA6*Y51)';
        pP46CPtildestart(:,ll) = (ones(1,c*n)*Y127*bA6*pbP6P4start(:,:,ll)*Y51)';
        pP46Chabitatstart(:,ll) = (ones(1,c*n)*Y127*pA6P4start*Y51)';

%         pP56CP(:,ll) = (ones(1,c*n)*Y127*bA6*Y51)';
        pP56CPtildestart(:,ll) = (ones(1,c*n)*Y127*bA6*pbP6P5start(:,:,ll)*Y51)';
        pP56Chabitatstart(:,ll) = (ones(1,c*n)*Y127*pA6P5start*Y51)';
        
        for ii = 1:c*n 
                % all contribution metrics will be nonzero in same place 
                    % corresponding to the states populated at the start of the
                    % annual cycle
                if pP46CPtildestart(ii,ll) == 0
                    P46elasCPtildestart(ii,pp,ll) = 0;
                else
                    P46elasCPtildestart(ii,pp,ll) = (d/pP46CPtildestart(ii,ll))*P46sensCPtildestart(ii,pp,ll);
                end

                if pP46Chabitatstart(ii,ll) == 0
                    P46elasChabitatstart(ii,pp,ll) = 0;
                else
                    P46elasChabitatstart(ii,pp,ll) = (d/pP46Chabitatstart(ii,ll))*P46sensChabitatstart(ii,pp,ll);
                end
                
                if pP56CPtildestart(ii,ll) == 0
                    P56elasCPtildestart(ii,pp,ll) = 0;
                else
                    P56elasCPtildestart(ii,pp,ll) = (d/pP56CPtildestart(ii,ll))*P56sensCPtildestart(ii,pp,ll);
                end

                if pP56Chabitatstart(ii,ll) == 0
                    P56elasChabitatstart(ii,pp,ll) = 0;
                else
                    P56elasChabitatstart(ii,pp,ll) = (d/pP56Chabitatstart(ii,ll))*P56sensChabitatstart(ii,pp,ll);
                end
        end
%         % Substituting all \delta for 0.01 in elasticities function
% %          P46SUBSelasCP = subs(P46elasCP,d,0.01);
%          P46SUBSelasCPtildestart = subs(P46elasCPtildestart,d,0.01);
%          P46SUBSelasChabitatstart = subs(P46elasChabitatstart,d,0.01);  
%          %
% %          P56SUBSelasCP = subs(P56elasCP,d,0.01);
%          P56SUBSelasCPtildestart = subs(P56elasCPtildestart,d,0.01);
%          P56SUBSelasChabitatstart = subs(P56elasChabitatstart,d,0.01);  

         % Substituting each \delta to be 1% of the life rate it is
         % being positively perturbed
         P46SUBSelasCPtildestart = subs(P46elasCPtildestart,d,deltasub);
         P46SUBSelasChabitatstart = subs(P46elasChabitatstart,d,deltasub); 
         %
         P56SUBSelasCPtildestart = subs(P56elasCPtildestart,d,deltasub);
         P56SUBSelasChabitatstart = subs(P56elasChabitatstart,d,deltasub); 
    end
end 

%% SENSITIVITY OF CONTRIBUTION METRICS WITH RESPECT TO PROPORTION MIGRATING V3
% MIGRATION HAPPENS BEFORE DEMOGRAPHY 
% there are 4 proportion rates for both stage 4 & 5, but within each stage,
% the proportion rates can only be perturbed in pairs, as the columns of
% P^i_k need to sum to 1. Thus there are only two available perturbations
% for each stage, and so 4 perturbations in total. 

% the sensitivity of \bC(\sP) is 0 for all as there is a specified path in
% season 6 (Z_{2,6}=0). season 6 is also when we are performing the
% perturbations. 


pmP46 = sym(zeros(c*n,c*n,height(PATH))); 
pmP56 = pmP46; 
% bold P terms are affected by perturbation to P^i_k matrices
pbP6P4 = pmP46; % path specified 
pbP6P5 = pbP6P4; % path specified
% formula for the contribution metrics following perturbation to P^4_6
pP46CP = sym(zeros(c*n,height(PATH)));
pP46CPtilde = pP46CP;
pP46Chabitat = pP46CP;
% formula for the contribution metrics following perturbation to P^5_6
pP56CP = sym(zeros(c*n,height(PATH)));
pP56CPtilde = pP56CP;
pP56Chabitat = pP56CP;

% sensitivity formula for perturbations to P^4_6
P46sensCPtilde = sym(zeros(c*n,2,height(PATH))); 
P46sensChabitat = P46sensCPtilde;
% elasticity formula for perturbations to P^4_6 life rates
P46elasCPtilde = P46sensCPtilde;
P46elasChabitat = P46elasCPtilde;
% sensitivity formula for perturbations to P^5_6 
P56sensCPtilde = sym(zeros(c*n,2,height(PATH))); 
P56sensChabitat = P56sensCPtilde;
% elasticity formula for perturbations to P^5_6 life rates
P56elasCPtilde = P56sensCPtilde;
P56elasChabitat = P56elasCPtilde;
% dvec(A_k)/dp^T is zero as does not contain P matrices 

for ll = 1:height(PATH)
    % constants that depend on which path is specified   
    EPP = zeros(n); EPP(PATH(ll,7),PATH(ll,6))=1; 
    bA6 = Ahat(:,:,6).*kron(EPP,ones(c));
    
    T_16 = zeros(c,c*n); T_16(:,(PATH(ll,7)-1)*c+1:PATH(ll,7)*c) = eye(c);
    T_26 = zeros(c*n,c); T_26((PATH(ll,6)-1)*c+1:PATH(ll,6)*c,:) = eye(c);
    
    for pp = 1:4
        syms d real % act as delta
        if pp == 1 
        % +ve perturbation is to P^4_6(3,2) and 
        % -ve perturbation is to P^4_6(4,2) similarly stage 5 equivalents
            ploc = [0 0 0 0; 0 0 0 0; 0 d 0 0; 0 -d 0 0];
            % 1% of life rate being positively perturbed for elasticity substitutions
            deltasub = 0.556*0.01;
        elseif pp == 2
        % +ve perturbation is to P^4_6(4,2) and 
        % -ve perturbation is to P^4_6(3,2) similarly stage 5 equivalents
            ploc = [0 0 0 0; 0 0 0 0; 0 -d 0 0; 0 d 0 0];
            % 1% of life rate being positively perturbed for elasticity substitutions
            deltasub = 0.444*0.01;
        elseif pp == 3
        % +ve perturbation is to P^4_6(3,3) and 
        % -ve perturbation is to P^4_6(4,3) similarly stage 5 equivalents
            ploc = [0 0 0 0; 0 0 0 0; 0 0 d 0; 0 0 -d 0];
            % 1% of life rate being positively perturbed for elasticity substitutions
            deltasub = 0.532*0.01;
        elseif pp==4 
        % +ve perturbation is to P^4_6(4,3) and 
        % -ve perturbation is to P^4_6(3,3) similarly stage 5 equivalents
            ploc = [0 0 0 0; 0 0 0 0; 0 0 -d 0; 0 0 d 0];
            % 1% of life rate being positively perturbed for elasticity substitutions
            deltasub = 0.468*0.01;
        end
                
        pP46 = P(:,3*n+1:4*n,6) + ploc; % perturbed version of P^4_6
        vecPp46 = reshape(pP46,n*n,1);
        
        pP56 = P(:,4*n+1:5*n,6) + ploc; % perturbed version of P^5_6
        vecPp56 = reshape(pP56,n*n,1);
        
        % perturbed mP_6
        for ii = 1:c 
            E = zeros(c); E(ii,ii) = 1;  
            if ii == 4 
                pmP46(:,:,ll) = pmP46(:,:,ll) + kron(pP46,E);
                pmP56(:,:,ll) = pmP56(:,:,ll) + kron(P(:,(ii-1)*n+1:ii*n,6),E);
            elseif ii == 5 
                pmP46(:,:,ll) = pmP46(:,:,ll) + kron(P(:,(ii-1)*n+1:ii*n,6),E);
                pmP56(:,:,ll) = pmP56(:,:,ll) + kron(pP56,E);
            else
                pmP46(:,:,ll) = pmP46(:,:,ll) + kron(P(:,(ii-1)*n+1:ii*n,6),E);
                pmP56(:,:,ll) = pmP56(:,:,ll) + kron(P(:,(ii-1)*n+1:ii*n,6),E);
            end
        end
        
        P46sensCPtilde(:,pp,ll) = X4*kronUU*kron(eye(c*n),bA6)*(Z46*X1*...
            kron(vecIn,(kron(T_26',T_16)*(X1*kron(diff(vecPp46),vecEc44)))));
        P46sensChabitat(:,pp,ll) = X4*kronYY*X2k(:,:,6)*(vecmS6.*(X1*kron(diff(vecPp46),vecEc44)));
        
        P56sensCPtilde(:,pp,ll) = X4*kronUU*kron(eye(c*n),bA6)*(Z46*X1*...
            kron(vecIn,(kron(T_26',T_16)*(X1*kron(diff(vecPp56),vecEc55)))));
        P56sensChabitat(:,pp,ll) = X4*kronYY*X2k(:,:,6)*(vecmS6.*(X1*kron(diff(vecPp56),vecEc55)));
        
        % perturbed contribution metrics
        pbP6P4(:,:,ll) = kron(eye(n),T_16*pmP46(:,:,ll)*T_26);
        pbP6P5(:,:,ll) = kron(eye(n),T_16*pmP56(:,:,ll)*T_26);
        
        pA6P4 = mD(:,:,6)*(mS(:,:,6).*pmP46(:,:,ll));
        pA6P5 = mD(:,:,6)*(mS(:,:,6).*pmP56(:,:,ll));
        
        pP46CP(:,ll) = (ones(1,c*n)*Y127*bA6*Y51)';
        pP46CPtilde(:,ll) = (ones(1,c*n)*Y127*bA6*pbP6P4(:,:,ll)*Y51)';
        pP46Chabitat(:,ll) = (ones(1,c*n)*Y127*pA6P4*Y51)';

        pP56CP(:,ll) = (ones(1,c*n)*Y127*bA6*Y51)';
        pP56CPtilde(:,ll) = (ones(1,c*n)*Y127*bA6*pbP6P5(:,:,ll)*Y51)';
        pP56Chabitat(:,ll) = (ones(1,c*n)*Y127*pA6P5*Y51)';
        
        for ii = 1:c*n 
                % all contribution metrics will be nonzero in same place 
                    % corresponding to the states populated at the start of the
                    % annual cycle
                if pP46CPtilde(ii,ll) == 0
                    P46elasCPtilde(ii,pp,ll) = 0;
                else
                    P46elasCPtilde(ii,pp,ll) = (d/pP46CPtilde(ii,ll))*P46sensCPtilde(ii,pp,ll);
                end

                if pP46Chabitat(ii,ll) == 0
                    P46elasChabitat(ii,pp,ll) = 0;
                else
                    P46elasChabitat(ii,pp,ll) = (d/pP46Chabitat(ii,ll))*P46sensChabitat(ii,pp,ll);
                end
                
                if pP56CPtilde(ii,ll) == 0
                    P56elasCPtilde(ii,pp,ll) = 0;
                else
                    P56elasCPtilde(ii,pp,ll) = (d/pP56CPtilde(ii,ll))*P56sensCPtilde(ii,pp,ll);
                end

                if pP56Chabitat(ii,ll) == 0
                    P56elasChabitat(ii,pp,ll) = 0;
                else
                    P56elasChabitat(ii,pp,ll) = (d/pP56Chabitat(ii,ll))*P56sensChabitat(ii,pp,ll);
                end
        end
%         % Substituting all \delta for 0.01 in elasticities function
%          P46SUBSelasCP = subs(P46elasCP,d,0.01);
%          P46SUBSelasCPtilde = subs(P46elasCPtilde,d,0.01);
%          P46SUBSelasChabitat = subs(P46elasChabitat,d,0.01);  
%          %
%          P56SUBSelasCP = subs(P56elasCP,d,0.01);
%          P56SUBSelasCPtilde = subs(P56elasCPtilde,d,0.01);
%          P56SUBSelasChabitat = subs(P56elasChabitat,d,0.01);  
         
         % Substituting each \delta to be 1% of the life rate it is
         % being positively perturbed
%          P46SUBSelasCP2 = subs(P46elasCP,d,deltasub);
         P46SUBSelasCPtilde2 = subs(P46elasCPtilde,d,deltasub);
         P46SUBSelasChabitat2 = subs(P46elasChabitat,d,deltasub); 
         %
%          P56SUBSelasCP2 = subs(P56elasCP,d,deltasub);
         P56SUBSelasCPtilde2 = subs(P56elasCPtilde,d,deltasub);
         P56SUBSelasChabitat2 = subs(P56elasChabitat,d,deltasub); 
    end
end 


%% APPLYING MULTIPLE PERTURBATIONS AT ONCE 

% perturbations to Pathway 2.
multpPATH = PATH(2,:); 
% constants that depend on which path is specified
EPP = zeros(n); EPP(multPATH(7),multPATH(6))=1; 
Z36 = reshape(kron(EPP,ones(c)),n*c*n*c,1);%constant Z_{3,k}for path ll
% used in calculation of dCPtilde/dp & depends which path is specified
bP6 = kron(eye(n),mP((multPATH(7)-1)*c+1:multPATH(7)*c,(multPATH(6)-1)*c+1:multPATH(6)*c,6));

% perturbation to proportion migrating rates (P^4_6 and P^5_6 affected.) 


% formula for the contribution metrics following all perturbations
pDSPCP = sym(zeros(c*n,1));
pDSPCPtilde = pDSPCP;
pDSPChabitat = pDSPCP;

% sensitivity formula for perturbations to D_{4,6}, S^4_6, S^5_6, P^4_6 and P^5_6
DSPsensCP = zeros(c*n,1); 
DSPsensCPtilde = DSPsensCP;
DSPsensChabitat = DSPsensCP;
% elasticity formula for perturbations to D_{4,6}, S^4_6, S^5_6, P^4_6 and P^5_6
DSPelasCP = sym(zeros(c*n,1)); 
DSPelasCPtilde = DSPelasCP;
DSPelasChabitat = DSPelasCP;

% constants used in all perturbations to D_46
E44 = zeros(n); E44(4,4)=1; % matrix E_{n,jj} j=4
vecE44 = reshape(E44,n*n,1);


% Carrying out perturbations
syms dD dS dP real % act as delta perturbations to D_{4,6}, S^{4 and 5}_6 and P^{4 and 5}_6

% perturbations to D matrix (D_46)
pmD6 = sym(zeros(c*n,c*n)); % only 1 path => no 3rd dimension
%
plocD = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 dD 0];
pD = D(:,(4-1)*c+1:4*c,6) + plocD; %perturbed version of D_46
vecDp = reshape(pD,c*c,1);
% perturbed mD_6
for jj = 1:n 
    if jj == 4 
        pmD6((jj-1)*c+1:jj*c,(jj-1)*c+1:jj*c,ll) = pD;
    else
        pmD6((jj-1)*c+1:jj*c,(jj-1)*c+1:jj*c,ll) = mD((jj-1)*c+1:jj*c,(jj-1)*c+1:jj*c,6);
    end
end

% perturbations to S matrices (S^4_6 and S^5_6)
pmS = sym(zeros(c*n,c*n)); 
%
plocS = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 dS 0 0];
pS46 = S(:,3*n+1:4*n,6) + plocS; % perturbed version of S^4_6
vecpS46 = reshape(pS46,n*n,1); 
pS56 = S(:,4*n+1:5*n,6) + plocS; % perturbed version of S^5_6
vecpS56 = reshape(pS56,n*n,1);
% perturbed mS_6
for ii = 1:c 
    E = zeros(c); E(ii,ii) = 1;  
    if ii == 4  % stage 4 perturbed
        pmS(:,:,ll) = pmS(:,:,ll) + kron(pS46,E);
    elseif ii == 5 % stage 5 perturbed
        pmS(:,:,ll) = pmS(:,:,ll) + kron(pS56,E);
    else
        pmS(:,:,ll) = pmS(:,:,ll) + kron(S(:,(ii-1)*n+1:ii*n,6),E);
    end
end
% parts within the sum over i. X1 is constant so can be pulled out
kronpS4 = kron(diff(vecpS46),vecEc44);
kronpS5 = kron(diff(vecpS56),vecEc55);

% perturbations to P matrices 
pmP = sym(zeros(c*n,c*n)); 
%
plocP = [0 0 0 0; 0 0 0 0; 0 dP 0 0; 0 -dP 0 0];
pP46 = P(:,3*n+1:4*n,6) + plocP; % perturbed version of P^4_6
vecpP46 = reshape(pP46,n*n,1); 
pP56 = P(:,4*n+1:5*n,6) + plocP; % perturbed version of P^5_6
vecpP56 = reshape(pP56,n*n,1);
% perturbed mP_6
for ii = 1:c 
    E = zeros(c); E(ii,ii) = 1;  
    if ii == 4  % stage 4 perturbed
        pmP(:,:,ll) = pmP(:,:,ll) + kron(pP46,E);
    elseif ii == 5 % stage 5 perturbed
        pmP(:,:,ll) = pmP(:,:,ll) + kron(pP56,E);
    else
        pmP(:,:,ll) = pmP(:,:,ll) + kron(P(:,(ii-1)*n+1:ii*n,6),E);
    end
end

dvecA6S = (X2k(:,:,6)*((X1*(kronpS4+kronpS5)).*vecmP6)).*Z36;
dvecA6Shabitat = X2k(:,:,6)*(X1*(kronpS4+kronpS5).*vecmP6);

% bold P terms are affected by perturbation to P^i_k matrices
pbP6P = pmP; % path specified 

SsensCP(:,pp,ll) = X4*kronYY*dvecA6S;
SsensCPtilde(:,pp,ll) = X4*kronUU*kron(bP6',eye(c*n))*dvecA6S;
SsensChabitat(:,pp,ll) = X4*kronYY*dvecA6Shabitat;

% perturbed contribution metrics
pbA6S(:,:,ll) = (mD(:,:,6)*pmS(:,:,ll)).*kron(EPP,ones(c)); % path specified
pA6S(:,:,ll) = mD(:,:,6)*(pmS(:,:,ll).*mP(:,:,6)); % no path specified

pSCP(:,ll) = (ones(1,c*n)*Y127*pbA6S(:,:,ll)*Y51)';
pSCPtilde(:,ll) = (ones(1,c*n)*Y127*pbA6S(:,:,ll)*bP6*Y51)';
pSChabitat(:,ll) = (ones(1,c*n)*Y127*pA6S(:,:,ll)*Y51)';

for ii = 1:c*n 
    % all contribution metrics will be nonzero in same place 
        % corresponding to the states populated at the start of the
        % annual cycle
    if pSCP(ii,ll) == 0 
        SelasCP(ii,pp,ll) = 0; % or NaN
    else
        SelasCP(ii,pp,ll) = (d/pSCP(ii,ll))*SsensCP(ii,pp,ll);
    end

    if pSCPtilde(ii,ll) == 0
        SelasCPtilde(ii,pp,ll) = 0;
    else
        SelasCPtilde(ii,pp,ll) = (d/pSCPtilde(ii,ll))*SsensCPtilde(ii,pp,ll);
    end

    if pSChabitat(ii,ll) == 0
        SelasChabitat(ii,pp,ll) = 0;
    else
        SelasChabitat(ii,pp,ll) = (d/pSChabitat(ii,ll))*SsensChabitat(ii,pp,ll);
    end
end

% Substituting all \delta for 0.01 in elasticities function
SSUBSelasCP = subs(SelasCP,d,0.01);
SSUBSelasCPtilde = subs(SelasCPtilde,d,0.01);
SSUBSelasChabitat = subs(SelasChabitat,d,0.01);  





%% HOW PERTURBATIONS TO DEMOGRAPHY AFFECT ASYMPTOTICS
% annual cycle matrix following perturbations to D_{3,6} life rates
% D36sA = zeros(c*n,c*n,7); % no path specified for \sA
D_36 = D(:,(3-1)*c+1:3*c,6);
D_46 = D(:,(4-1)*c+1:4*c,6);
D36mD = mD(:,:,6); % to store perturbed version of \mD_6
D46mD = mD(:,:,6);
lambdaD36sA = zeros(100,7);
lambdaD46sA = zeros(100,7);
senslambdaD36 = zeros(1,7);
senslambdaD46 = zeros(1,7);
elaslambdaD36 = zeros(1,7);
elaslambdaD46 = zeros(1,7);
for pp = 1:7  %number of non-zero entries in D_36 and D_46
        if pp==1 % perturbation is to D_j6(1,3)
            d = linspace(-2*268,2*268);% possible values of delta_1
            for dd = 1:100 
                pD_36 = D_36; 
                pD_36(1,3) = D_36(1,3) + d(dd); %perturbed version of D_36
                D36mD((3-1)*c+1:3*c,(3-1)*c+1:3*c) = pD_36;
                D36sA = Y127*D36mD*mM(:,:,6)*Y51;
                lambdaD36sA(dd,pp) = max(abs(eig(D36sA)));
                
                pD_46 = D_46; 
                pD_46(1,3) = D_46(1,3) + d(dd); %perturbed version of D_46
                D46mD((4-1)*c+1:4*c,(4-1)*c+1:4*c) = pD_46;
                D46sA = Y127*D46mD*mM(:,:,6)*Y51;
                lambdaD46sA(dd,pp) = max(abs(eig(D46sA)));
            end
        elseif pp==2 % perturbation is to D_j6(1,4)
            d = linspace(-2*268,2*268);% possible values of delta_2
            for dd = 1:100 
                pD_36 = D_36; 
                pD_36(1,4) = D_36(1,4) + d(dd); %perturbed version of D_36
                D36mD((3-1)*c+1:3*c,(3-1)*c+1:3*c) = pD_36;
                D36sA = Y127*D36mD*mM(:,:,6)*Y51;
                lambdaD36sA(dd,pp) = max(abs(eig(D36sA)));
                
                pD_46 = D_46; 
                pD_46(1,4) = D_46(1,4) + d(dd); %perturbed version of D_46
                D46mD((4-1)*c+1:4*c,(4-1)*c+1:4*c) = pD_46;
                D46sA = Y127*D46mD*mM(:,:,6)*Y51;
                lambdaD46sA(dd,pp) = max(abs(eig(D46sA)));
            end
        elseif pp==3    % perturbation is to D_j6(1,5)
            d = linspace(-2*89,2*89);% possible values of delta_3
            for dd = 1:100
                pD_36 = D_36; 
                pD_36(1,5) = D_36(1,5) + d(dd); %perturbed version of D_36
                D36mD((3-1)*c+1:3*c,(3-1)*c+1:3*c) = pD_36;
                D36sA = Y127*D36mD*mM(:,:,6)*Y51;
                lambdaD36sA(dd,pp) = max(abs(eig(D36sA)));
                
                pD_46 = D_46; 
                pD_46(1,5) = D_46(1,5) + d(dd); %perturbed version of D_46
                D46mD((4-1)*c+1:4*c,(4-1)*c+1:4*c) = pD_46;
                D46sA = Y127*D46mD*mM(:,:,6)*Y51;
                lambdaD46sA(dd,pp) = max(abs(eig(D46sA)));
            end
        elseif pp==4    % perturbation is to D_j6(4,1)
            d = linspace(-0.416,0.584);% possible values of delta_4
            for dd = 1:100 
                pD_36 = D_36; 
                pD_36(4,1) = D_36(4,1) + d(dd); %perturbed version of D_36
                D36mD((3-1)*c+1:3*c,(3-1)*c+1:3*c) = pD_36;
                D36sA = Y127*D36mD*mM(:,:,6)*Y51;
                lambdaD36sA(dd,pp) = max(abs(eig(D36sA)));
                
                pD_46 = D_46; 
                pD_46(4,1) = D_46(4,1) + d(dd); %perturbed version of D_46
                D46mD((4-1)*c+1:4*c,(4-1)*c+1:4*c) = pD_46;
                D46sA = Y127*D46mD*mM(:,:,6)*Y51;
                lambdaD46sA(dd,pp) = max(abs(eig(D46sA)));
            end
        elseif pp==5    % perturbation is to D_j6(5,3)
            d = linspace(-0.56,0.44);% possible values of delta_5
            for dd = 1:100 
                pD_36 = D_36; 
                pD_36(5,3) = D_36(5,3) + d(dd); %perturbed version of D_36
                D36mD((3-1)*c+1:3*c,(3-1)*c+1:3*c) = pD_36;
                D36sA = Y127*D36mD*mM(:,:,6)*Y51;
                lambdaD36sA(dd,pp) = max(abs(eig(D36sA)));
                
                pD_46 = D_46; 
                pD_46(5,3) = D_46(5,3) + d(dd); %perturbed version of D_46
                D46mD((4-1)*c+1:4*c,(4-1)*c+1:4*c) = pD_46;
                D46sA = Y127*D46mD*mM(:,:,6)*Y51;
                lambdaD46sA(dd,pp) = max(abs(eig(D46sA)));
            end
        elseif pp==6    % perturbation is to D_j6(5,4)
            d = linspace(-0.56,0.44);% possible values of delta_6
            for dd = 1:100 
                pD_36 = D_36; 
                pD_36(5,4) = D_36(5,4) + d(dd); %perturbed version of D_36
                D36mD((3-1)*c+1:3*c,(3-1)*c+1:3*c) = pD_36;
                D36sA = Y127*D36mD*mM(:,:,6)*Y51;
                lambdaD36sA(dd,pp) = max(abs(eig(D36sA)));
                
                pD_46 = D_46; 
                pD_46(5,4) = D_46(5,4) + d(dd); %perturbed version of D_46
                D46mD((4-1)*c+1:4*c,(4-1)*c+1:4*c) = pD_46;
                D46sA = Y127*D46mD*mM(:,:,6)*Y51;
                lambdaD46sA(dd,pp) = max(abs(eig(D46sA)));
            end
        elseif pp==7    % perturbation is to D_j6(5,5)
            d = linspace(-0.17,0.83);% possible values of delta_4
            for dd = 1:100 
                pD_36 = D_36; 
                pD_36(5,5) = D_36(5,5) + d(dd); %perturbed version of D_36
                D36mD((3-1)*c+1:3*c,(3-1)*c+1:3*c) = pD_36;
                D36sA = Y127*D36mD*mM(:,:,6)*Y51;
                lambdaD36sA(dd,pp) = max(abs(eig(D36sA)));
                
                pD_46 = D_46; 
                pD_46(5,5) = D_46(5,5) + d(dd); %perturbed version of D_46
                D46mD((4-1)*c+1:4*c,(4-1)*c+1:4*c) = pD_46;
                D46sA = Y127*D46mD*mM(:,:,6)*Y51;
                lambdaD46sA(dd,pp) = max(abs(eig(D46sA)));
            end
        end
        
        if pp == 1 || pp == 2 || pp == 3    % reproductive life rates 
            figure(1)
            subplot(1,2,1)
            plot(d,lambdaD36sA(:,pp),'LineWidth',2)
            hold on
            plot(d,lambdaD46sA(:,pp),'--','LineWidth',2)
        else 
            figure(1)
            subplot(1,2,2)
            plot(d,lambdaD36sA(:,pp),'LineWidth',2)
            hold on
            plot(d,lambdaD46sA(:,pp),'--','LineWidth',2)
        end 
        
        % find entries of d that are closest to 0 (both +ve and -ve)
        nve0 = max(d(d<=0)); 
        nve0 = find(abs(d-nve0)<0.00001);
        pve0 = min(d(d>=0));
        pve0 = find(abs(d-pve0)<0.00001);
        % estimate of slope at delta = 0 
        senslambdaD36(pp) = (lambdaD36sA(max(nve0,pve0),pp)-lambdaD36sA(min(nve0,pve0),pp))/(d(max(nve0,pve0))-d(min(nve0,pve0)));
        senslambdaD46(pp) = (lambdaD46sA(max(nve0,pve0),pp)-lambdaD46sA(min(nve0,pve0),pp))/(d(max(nve0,pve0))-d(min(nve0,pve0)));

%         elaslambdaD36(pp) = (d/Cptilde(ll,ii))*P56sensCPtildestart(ii,pp,ll);        
end

figure(1)
subplot(1,2,1)
legend('delta1D36','delta1D46','delta2D36','delta2D46','delta3D36','delta3D46')
hold on
subplot(1,2,2)
legend('delta4D36','delta4D46','delta5D36','delta5D46','delta6D36','delta6D46','delta7D36','delta7D46')

%% HOW PERTRUBATIONS TO MOVEMENT SURVIVAL AFFECT ASYMPTOTICS 
% annual cycle matrix following perturbations to S^i_6, where i \in {4,5}
S_46 = S(:,(4-1)*n+1:4*n,6);
S_56 = S(:,(5-1)*n+1:5*n,6);
S46mS = mS(:,:,6); % to store perturbed version of \mS_6
S56mS = mS(:,:,6);
lambdaS46sA = zeros(100,3); % there are 3 possible entries to perturb
lambdaS56sA = zeros(100,3);
senslambdaS46 = zeros(1,3);
senslambdaS56 = zeros(1,3);

for pp = 1:3  %number of entries in S^4_6 and S^5_6 subject to perturbation
        if pp==1 % perturbation is to S^i_6(3,2)
            d = linspace(-0.733,0.267);% possible values of delta_1
            for dd = 1:100 
                % perturb the (i,i)-th entry of (3,2)_th block of mS_6. i=4
                pS46mS = S46mS;
                pS46mS((3-1)*c+4,(2-1)*c+4) = S46mS((3-1)*c+4,(2-1)*c+4) + d(dd); 
                S46sA = Y127*mD(:,:,6)*(pS46mS.*mP(:,:,6))*Y51;
                lambdaS46sA(dd,pp) = max(abs(eig(S46sA)));
                
                % perturb the (i,i)-th entry of (3,2)_th block of mS_6. i=5
                pS56mS = S56mS;
                pS56mS((3-1)*c+5,(2-1)*c+5) = S56mS((3-1)*c+5,(2-1)*c+5) + d(dd); 
                S56sA = Y127*mD(:,:,6)*(pS56mS.*mP(:,:,6))*Y51;
                lambdaS56sA(dd,pp) = max(abs(eig(S56sA)));
            end
        elseif pp==2 % perturbation is to S^i_6(4,2)
            d = linspace(-0.544,0.456);% possible values of delta_2
            for dd = 1:100
                pS46mS = S46mS;
                % perturb the (i,i)-th entry of (4,2)_th block of mS_6. i=4
                pS46mS((4-1)*c+4,(2-1)*c+4) = S46mS((4-1)*c+4,(2-1)*c+4) + d(dd); 
                S46sA = Y127*mD(:,:,6)*(pS46mS.*mP(:,:,6))*Y51;
                lambdaS46sA(dd,pp) = max(abs(eig(S46sA)));
                
                % perturb the (i,i)-th entry of (4,2)_th block of mS_6. i=5
                pS56mS = S56mS;
                pS56mS((4-1)*c+5,(2-1)*c+5) = S56mS((4-1)*c+5,(2-1)*c+5) + d(dd); 
                S56sA = Y127*mD(:,:,6)*(pS56mS.*mP(:,:,6))*Y51;
                lambdaS56sA(dd,pp) = max(abs(eig(S56sA)));
            end
        elseif pp==3    % perturbation is to S^i_6(4,3)
            d = linspace(-0.742,0.258);% possible values of delta_3
            for dd = 1:100
                pS46mS = S46mS;
                % perturb the (i,i)-th entry of (4,3)_th block of mS_6. i=4
                pS46mS((4-1)*c+4,(3-1)*c+4) = S46mS((4-1)*c+4,(3-1)*c+4) + d(dd); 
                S46sA = Y127*mD(:,:,6)*(pS46mS.*mP(:,:,6))*Y51;
                lambdaS46sA(dd,pp) = max(abs(eig(S46sA)));
                
                % perturb the (i,i)-th entry of (4,3)_th block of mS_6. i=5
                pS56mS = S56mS;
                pS56mS((4-1)*c+5,(3-1)*c+5) = S56mS((4-1)*c+5,(3-1)*c+5) + d(dd); 
                S56sA = Y127*mD(:,:,6)*(pS56mS.*mP(:,:,6))*Y51;
                lambdaS56sA(dd,pp) = max(abs(eig(S56sA)));
            end
        end
        
        figure(2)
        plot(d,lambdaS46sA(:,pp),'LineWidth',2)
        hold on
        plot(d,lambdaS56sA(:,pp),'--','LineWidth',2)
        hold on
        
        % find the negative and positive entries of d that are closest to 0
        nve0 = max(d(d<=0)); 
        nve0 = find(abs(d-nve0)<0.00001);
        pve0 = min(d(d>=0));
        pve0 = find(abs(d-pve0)<0.00001);
        % estimate of slope at delta = 0 
        senslambdaS46(pp) = (lambdaS46sA(max(nve0,pve0),pp)-lambdaS46sA(min(nve0,pve0),pp))/(d(max(nve0,pve0))-d(min(nve0,pve0)));
        senslambdaS56(pp) = (lambdaS56sA(max(nve0,pve0),pp)-lambdaS56sA(min(nve0,pve0),pp))/(d(max(nve0,pve0))-d(min(nve0,pve0)));
        
end

figure(2)
legend('delta1S46','delta1S56','delta2S46','delta2S56','delta3S46','delta3S56')

%% HOW PERTRUBATIONS TO MOVEMENT PROPORTION AFFECT ASYMPTOTICS 
% annual cycle matrix following perturbations to P^i_6, where i \in {4,5}
P_46 = P(:,(4-1)*n+1:4*n,6);
P_56 = P(:,(5-1)*n+1:5*n,6);
P46mP = mP(:,:,6); % to store perturbed version of \mP_6
P56mP = mP(:,:,6);
lambdaP46sA = zeros(100,2); %there are 4 possible entries to perturb but
lambdaP56sA = zeros(100,2);  %two rates must be perturbed at once to ensure
senslambdaP46 = zeros(1,2);  %proportions sum to 1, so have 2 deltas
senslambdaP56 = zeros(1,2);

for pp = 1:2  %number of deltas that perturb P^4_6 and P^5_6
        if pp==1 % perturbation is to P^i_6(3,2) and P^i_^(4,2)
            d = linspace(-0.566,0.444);% possible values of delta_1
            for dd = 1:100 
                % +ve perturbation to the (i,i)-th entry of (3,2)_th block
                % of mP_6. i=4
                pP46mP = P46mP;
                pP46mP((3-1)*c+4,(2-1)*c+4) = P46mP((3-1)*c+4,(2-1)*c+4) + d(dd);
                % -ve perturbation to the (i,i)-th entry of (4,2)_th block
                % of mP_6. i=4
                pP46mP((4-1)*c+4,(2-1)*c+4) = P46mP((4-1)*c+4,(2-1)*c+4) - d(dd);
                P46sA = Y127*mD(:,:,6)*(mS(:,:,6).*pP46mP)*Y51;
                lambdaP46sA(dd,pp) = max(abs(eig(P46sA)));
                
                % +ve perturbation to the (i,i)-th entry of (3,2)_th block
                % of mP_6. i=5
                pP56mP = P56mP;
                pP56mP((3-1)*c+5,(2-1)*c+5) = P56mP((3-1)*c+5,(2-1)*c+5) + d(dd);
                % -ve perturbation to the (i,i)-th entry of (4,2)_th block
                % of mP_6. i=5
                pP56mP((4-1)*c+5,(2-1)*c+5) = P56mP((4-1)*c+5,(2-1)*c+5) - d(dd);
                P56sA = Y127*mD(:,:,6)*(mS(:,:,6).*pP56mP)*Y51;
                lambdaP56sA(dd,pp) = max(abs(eig(P56sA)));
            end
        elseif pp==2 % perturbation is to P^i_6(3,3) and P^i_^(4,3)
            d = linspace(-0.532,0.468);% possible values of delta_2
            for dd = 1:100 
                % +ve perturbation to the (i,i)-th entry of (3,3)_th block
                % of mP_6. i=4
                pP46mP = P46mP;
                pP46mP((3-1)*c+4,(3-1)*c+4) = P46mP((3-1)*c+4,(3-1)*c+4) + d(dd);
                % -ve perturbation to the (i,i)-th entry of (4,3)_th block
                % of mP_6. i=4
                pP46mP((4-1)*c+4,(3-1)*c+4) = P46mP((4-1)*c+4,(3-1)*c+4) - d(dd);
                P46sA = Y127*mD(:,:,6)*(mS(:,:,6).*pP46mP)*Y51;
                lambdaP46sA(dd,pp) = max(abs(eig(P46sA)));
                
                % +ve perturbation to the (i,i)-th entry of (3,3)_th block
                % of mP_6. i=5
                pP56mP = P56mP;
                pP56mP((3-1)*c+5,(3-1)*c+5) = P56mP((3-1)*c+5,(3-1)*c+5) + d(dd);
                % -ve perturbation to the (i,i)-th entry of (4,3)_th block
                % of mP_6. i=5
                pP56mP((4-1)*c+5,(3-1)*c+5) = P56mP((4-1)*c+5,(3-1)*c+5) - d(dd);
                P56sA = Y127*mD(:,:,6)*(mS(:,:,6).*pP56mP)*Y51;
                lambdaP56sA(dd,pp) = max(abs(eig(P56sA)));
            end
        end
        
        figure(3)
        plot(d,lambdaP46sA(:,pp),'LineWidth',2)
        hold on
        plot(d,lambdaP56sA(:,pp),'--','LineWidth',2)
        hold on
        
        % find the negative and positive entries of d that are closest to 0
        nve0 = max(d(d<=0)); 
        nve0 = find(abs(d-nve0)<0.00001);
        pve0 = min(d(d>=0));
        pve0 = find(abs(d-pve0)<0.00001);
        % estimate of slope at delta = 0 
        senslambdaP46(pp) = (lambdaP46sA(max(nve0,pve0),pp)-lambdaP46sA(min(nve0,pve0),pp))/(d(max(nve0,pve0))-d(min(nve0,pve0)));
        senslambdaP56(pp) = (lambdaP56sA(max(nve0,pve0),pp)-lambdaP56sA(min(nve0,pve0),pp))/(d(max(nve0,pve0))-d(min(nve0,pve0)));
        
end

figure(3)
legend('delta1P46','delta1P56','deltaPS46','deltaPS56')


%%
% %% PLOTS FOR PAPER
%  
% % Perturbations to demography survival rates (delta_4,...,delta_7) 
% figure(1)
% for ll = 1:4
%     subplot(2,2,ll);  % each subplot is associated with a different pathway
%     if ll == 1 || ll == 3   % demography happens in habitat 3
%         DsensCP = D36sensCP;
%         DsensCPtilde = D36sensCPtilde;
%         DsensChabitat = D36sensChabitat;
%     else                    % demography happens in habitat 4
%         DsensCP = D46sensCP;
%         DsensCPtilde = D46sensCPtilde;
%         DsensChabitat = D46sensChabitat;
%     end 
%     for pp = 4:7    % delta for perturbations to demography survival rates
%         if pp == 4
%             % x is biologically reasonable values
%             x = linspace(-0.416,0.584); 
%             y = Cp(ll,2) + DsensCP(2,pp,ll)*x;
%             plot(x,y,"r-",'LineWidth',2);
%             hold on
%             y = Cptilde(ll,2) + DsensCPtilde(2,pp,ll)*x;
%             plot(x,y,"r--",'LineWidth',2);
%             hold on
%             y = ChPtilde(1,2) + DsensChabitat(2,pp,ll)*x;
%             plot(x,y,"r:",'LineWidth',2);
%             hold on
%         elseif pp == 5 
%             % x is biologically reasonable values
%             x = linspace(-0.56,0.44); 
%             y = Cp(ll,2) + DsensCP(2,pp,ll)*x;
%             plot(x,y,"b-",'LineWidth',2);
%             hold on
%             y = Cptilde(ll,2) + DsensCPtilde(2,pp,ll)*x;
%             plot(x,y,"b--",'LineWidth',2);
%             hold on
%             y = ChPtilde(1,2) + DsensChabitat(2,pp,ll)*x;
%             plot(x,y,"b:",'LineWidth',2);
%             hold on
%         elseif pp == 6
%             % x is biologically reasonable values
%             x = linspace(-0.56,0.44); 
%             y = Cp(ll,2) + DsensCP(2,pp,ll)*x;
%             plot(x,y,"g-",'LineWidth',2);
%             hold on
%             y = Cptilde(ll,2) + DsensCPtilde(2,pp,ll)*x;
%             plot(x,y,"g--",'LineWidth',2);
%             hold on
%             y = ChPtilde(1,2) + DsensChabitat(2,pp,ll)*x;
%             plot(x,y,"g:",'LineWidth',2);
%             hold on
%         elseif pp == 7
%             % x is biologically reasonable values
%             x = linspace(-0.17,0.83); 
%             y = Cp(ll,2) + DsensCP(2,pp,ll)*x;
%             plot(x,y,"y-",'LineWidth',2);
%             hold on
%             y = Cptilde(ll,2) + DsensCPtilde(2,pp,ll)*x;
%             plot(x,y,"y--",'LineWidth',2);
%             hold on
%             y = ChPtilde(1,2) + DsensChabitat(2,pp,ll)*x;
%             plot(x,y,"y:",'LineWidth',2);
%             hold on
%         end
%     end
% end
% 
% % Perturbations to demography survival rates (delta_4,...,delta_7) 
% figure(2)
% for pp = 1:4    % count through delta_4...delta_7
%     subplot(2,2,pp);  % each subplot is associated with a different delta
%     % set x to be biologically reasonable values
%     if pp == 1
%         x = linspace(-0.416,0.584);
%     elseif pp == 2 || pp == 3
%         x = linspace(-0.56,0.44); 
%     elseif pp == 4
%         x = linspace(-0.17,0.83); 
%     end 
%     
%     for ll = 1:4 % count through different pathways
%         if ll == 1   % demography happens in habitat 3
%             y = Cp(ll,2) + D36sensCP(2,pp,ll)*x;
%             plot(x,y,'r-','LineWidth',2)
%             hold on
%             y = Cptilde(ll,2) + D36sensCPtilde(2,pp,ll)*x;
%             plot(x,y,"y-",'LineWidth',2);
%             hold on
%         elseif ll == 2  % demography happens in habitat 4
%             y = Cp(ll,2) + D46sensCP(2,pp,ll)*x;
%             plot(x,y,'r--','LineWidth',2);
%             hold on
%             y = Cptilde(ll,2) + D46sensCPtilde(2,pp,ll)*x;
%             plot(x,y,"y--",'LineWidth',2);
%             hold on
%         elseif ll == 3  % demography happens in habitat 3
%             y = Cp(ll,2) + D36sensCP(2,pp,ll)*x;
%             plot(x,y,'r:','LineWidth',2)
%             hold on
%             y = Cptilde(ll,2) + D36sensCPtilde(2,pp,ll)*x;
%             plot(x,y,"y:",'LineWidth',2);
%             hold on
%         elseif ll == 4  % demography happens in habitat 4
%             y = Cp(ll,2) + D46sensCP(2,pp,ll)*x;
%             plot(x,y,'r-.','LineWidth',2);
%             hold on
%             y = Cptilde(ll,2) + D46sensCPtilde(2,pp,ll)*x;
%             plot(x,y,"y-.",'LineWidth',2);
%             hold on
%         end 
%         y = ChPtilde(1,2) + DsensChabitat(2,pp,ll)*x;
%         plot(x,y,"g",'LineWidth',2);
%         hold on
%     end
% end
% 
%         
%         
% 
% 
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXTRA CODE RELATING TO MONARCH MODEL
% % SPECIFYING A PATH IN JUST ONE SEASON --- RUN THROUGH ALL OPTIONS
% Phi = 0; 
% PATH = zeros(28,s+1); %28 distict paths
% % populate PATH1 with all distinct paths
% PATH(1,1:2) = 1; % all valid season 1 paths
% PATH(2,2:3) = 1; % all valid season 2 paths
% PATH(3,3:4) = 1; % all valid season 3 paths
% PATH(4,4)=1; PATH(4,5)=2; % all valid season 4 paths
% PATH(5:6,5:6)=2; PATH(6,6)=3; % all valid season 5 paths
% PATH(7:8,6)=2; PATH(9:10,6)=3; PATH([7,9],7)=3; PATH([8,10],7)=4; % all valid season 6 paths
% PATH(11:12,7)=3; PATH(13:14,7)=4; PATH([11,13],8)=3; PATH([12,14],8)=4; % all valid season 7 paths
% PATH(15:16,8)=3; PATH(17:18,8)=4; PATH([15,17],9)=3; PATH([16,18],9)=4; % all valid season 8 paths
% PATH(19:20,9)=3; PATH(21:22,9)=4; PATH([19,21],10)=2; PATH([20,22],10)=3; % all valid season 9 paths
% PATH(23,10)=2; PATH(24:25,10)=3; PATH([23,25],11)=2; PATH(24,11)=1; % all valid season 10 paths
% PATH(26:27,11:12)=1; PATH(27,11)=2; % all valid season 11 paths
% PATH(28,12:13)=1; % all valid season 12 paths
% 
% 
% [row,col] = find(mD(:,:,kk));  %locate the non-zero entries of mD


%         if SS(PATH(ll,7),PATH(ll,6)) == 0 %if path does not exist
%             S16sensCP(:,ll) = NaN; %sensitivities are 0
%         else
%             syms d %act as delta
%             ploc = 
%             SS(PATH(ll,7),PATH(ll,6)) = SS(
%             S16sensCP(:,ll) = X4*kronYY*((X2(:,:,6)*



% % annual cycle matrix following perturbations to D_{3,6} life rates
% D46sA = zeros(c*n,c*n,7); % no path specified for \sA
% D46mD = mD(:,:,6); % to store perturbed version of \mD_6
% lambdaD46sA = zeros(1,7);
% for pp = 1:7 
%     %number of non-zero entries in D_46
%         syms d      %act as delta 
%         if pp==1 % perturbation is to D_46(1,3)
%             % location of perturbation
%             ploc = [0 0 d 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]; 
%         elseif pp==2 % perturbation is to D_46(1,4)
%             % location of perturbation
%             ploc = [0 0 0 d 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
%         elseif pp==3    % perturbation is to D_46(1,5)
%             % location of perturbation
%             ploc = [0 0 0 0 d; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
%         elseif pp==4    % perturbation is to D_46(4,1)
%             % location of perturbation
%             ploc = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; d 0 0 0 0; 0 0 0 0 0];
%         elseif pp==5    % perturbation is to D_46(5,3)
%             % location of perturbation
%             ploc = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 d 0 0];
%         elseif pp==6    % perturbation is to D_46(5,4)
%             % location of perturbation
%             ploc = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 d 0];
%         elseif pp==7    % perturbation is to D_46(5,5)
%             % location of perturbation
%             ploc = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 d];
%         end
%         
%         pD_46 = D_46 + ploc; %perturbed version of D_46
%         D46mD((4-1)*c+1:4*c,(4-1)*c+1:4*c,6) = pD_46;
%         
%         D46sA(:,:,pp) = Y127*D46mD*mM(:,:,6)*Y51;
%         
%         lambdaD46sA = max(abs(eig(D46sA(:,:,pp))));
% end


% % elasticities formula using diagonal matrices
% D36elasCP(:,pp,ll) = inv(diagCp(:,:,ll))*D36sensCP(:,pp,ll)*d;
% D36elasCPtilde(:,pp,ll) = inv(diagCptilde(:,:,ll))*D36sensCPtilde(:,pp,ll)*d;
% D36elasChabitat(:,pp,ll) = inv(diagChabitat(:,:,ll))*D36sensChabitat(:,pp,ll)*d;