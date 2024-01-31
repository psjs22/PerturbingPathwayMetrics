%%%%%%%%%%%%%%%%%%%%%%%%% SENSITIVITY OF LAMBDA  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% MONARCH BUTTERFLY EXAMPLE %%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
% Compute matrices for model mm in {1,2,3,4,5,6,7,8,9} 
mm = 6; % specify migratory model
[A,Ahat,n,c,s,mP,D,P,S,M,mD,mS,mM] = MigModel(mm);

AA = eye(c*n);
for kk = 1:s
    AA = A(:,:,kk)*AA;  % full annual cycle matrix, \sA
end 

% unperturbed lambda and corresponging eigenvectors
[W,lambda,V] = eig(AA);
maxpos = find(max(abs(lambda)));
V = V(:,maxpos);    % dominant left eigenvector 
lambda = lambda(maxpos,maxpos); % dominant eigenvalue
W = W(:,maxpos);    % dominant right eigenvector

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

X2k = zeros((c*n)^2,(c*n)^2,s); % constant X_{2,k}
for kk=1:s
    X2k(:,:,kk) = kron(eye(c*n),mD(:,:,kk));
end

Ec44 = zeros(c); Ec44(4,4)=1; %matrix E_{c,ii} i=4
vecEc44 = reshape(Ec44,c*c,1);

Ec55 = zeros(c); Ec55(5,5)=1; %matrix E_{c,ii} i=5
vecEc55 = reshape(Ec55,c*c,1);

% constant for path NOT specified in season 6
Z1k = zeros((c*n)^2,(c*n)^2,s); %constant Z_{1,k}
for kk = 1:s
    % all seasons do not have a specified path
    Z1k(:,:,kk) = kron(mM(:,:,kk)',eye(c*n));
end

Z2k = X2k; % no paths specified

Z3k = ones(c*n*c*n,1); % no paths specified

Y51 = eye(c*n); % to store constant Y^5_1
Y127 = Y51; % to store constant Y^12_7
for kk=1:5   
    Y51 = A(:,:,kk)*Y51;
    Y127 = A(:,:,kk+6)*Y127;    %only goes up to season 11
end
Y127 = A(:,:,12)*Y127; 
kronYY = kron(Y51',Y127); % (Y^5_1)^T \otimes Y^12_7 

vecmP6 = reshape(mP(:,:,6),c*n*c*n,1);
vecmS6 = reshape(mS(:,:,6),c*n*c*n,1);

vecIn = reshape(eye(n),n*n,1);

%% SENSITIVITY OF LAMBDA WITH RESPECT TO DEMOGRAPHIC RATES
% MIGRATION HAPPENS BEFORE DEMOGRAPHY => PATH(i,7) indicates which
% habitat's life rates will influence the Cmetrics

% for PATH([1,3],:) interested in D_{3,6}
D_36 = D(:,(3-1)*c+1:3*c,6);
% for PATH([2,4],:) interested in D_{4,6}
D_46 = D(:,(4-1)*c+1:4*c,6);

% sensitivity formula for perturbations to D_{3,6} life rates
D36sensLambda = sym(zeros(1,7)); % 7 is number of non zero life rates
% elasticity formula for perturbations to D_{3,6} life rates
D36elasLambda = sym(zeros(1,7)); % 7 is number of non zero life rates
% perturbing life rates in D_36 one at a time
% constants used in all perturbations to D_36
E33 = zeros(n); E33(3,3)=1; % matrix E_{n,jj} j=3
vecE33 = reshape(E33,n*n,1);
   
% compute sensitivities for all possible perturbations to D_36 
for pp = 1:7 %number of non-zero entries in D_36
    syms d      %act as delta 
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

    dvecA6D3 = (Z1k(:,:,6)*X1*kron(vecE33,diff(vecDp_36))).*Z3k;
    dvecAAD3 = kronYY*dvecA6D3;

    D36sensLambda(:,pp) = (kron(W',V')/(V'*W))*dvecAAD3;

    D36elasLambda(:,pp) = d/lambda*D36sensLambda(:,pp);
    D36SUBSelasLambda = subs(D36elasLambda,d,deltasub);
end            

%%
% sensitivity formula for pertrubations to D_{4,6} life rates 
D46sensLambda = sym(zeros(1,7)); % 7 is number of non zero life rates
% elasticity formula for perturbations to D_{3,6} life rates
D46elasLambda = sym(zeros(1,7)); % 7 is number of non zero life rates
% perturbing life rates in D_36 one at a time
% constants used in all perturbations to D_36
E44 = zeros(n); E44(4,4)=1; % matrix E_{n,jj} j=4
vecE44 = reshape(E44,n*n,1);
   
% compute sensitivities for all possible perturbations to D_46 
for pp = 1:7 %number of non-zero entries in D_46
    syms d      %act as delta 
    if pp==1 % perturbation is to D_46(1,3)
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

    dvecA6D4 = (Z1k(:,:,6)*X1*kron(vecE44,diff(vecDp_46))).*Z3k;
    dvecAAD4 = kronYY*dvecA6D4;

    D46sensLambda(:,pp) = (kron(W',V')/(V'*W))*dvecAAD4;

    D46elasLambda(:,pp) = d/lambda*D46sensLambda(:,pp);
    D46SUBSelasLambda = subs(D46elasLambda,d,deltasub);
end       

%% SENSITIVITY OF LAMBDA WITH RESPECT TO MOVEMENT SURVIVAL RATES
% MIGRATION HAPPENS BEFORE DEMOGRAPHY
% assume stages 1, 2 and 3 cannnot alter their movement survival rates
% there are 3 life rates that could change in both stage 4 and 5, so we
% will do a total of 6 perturbations

% sensitivity formula for perturbations to S^4_6 
S46sensLambda = sym(zeros(1,3)); 
% sensitivity formula for perturbations to S^5_6 
S56sensLambda = sym(zeros(1,3));
% elasticity formula for perturbations to S^4_6 life rates
S46elasLambda = sym(zeros(1,3));
% elasticity formula for perturbations to S^5_6 life rates
S56elasLambda = sym(zeros(1,3)); 
% perturbing movement rates associated with specified paths if non-zero
%%% assumption that new paths cannot be created 

% compute sensitivities for all possible perturbations to S^4_6 and S^5_6
for pp = 1:6
    syms d % act as delta
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

        dvecA6S4 = (X2k(:,:,6)*((X1*kron(diff(vecSp46),vecEcii)).*vecmP6)).*Z3k;
        dvecAAS4 = kronYY*dvecA6S4;

        S46sensLambda(:,pp) = (kron(W',V')/(V'*W))*dvecAAS4;
        S46elasLambda(:,pp) = d/lambda*S46sensLambda(:,pp);
        S46SUBSelasLambda = subs(S46elasLambda,d,deltasub);

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

        dvecA6S5 = (X2k(:,:,6)*((X1*kron(diff(vecSp56),vecEcii)).*vecmP6)).*Z3k;
        dvecAAS5 = kronYY*dvecA6S5;

        S56sensLambda(:,pp-3) = (kron(W',V')/(V'*W))*dvecAAS5;
        S56elasLambda(:,pp-3) = d/lambda*S56sensLambda(:,pp-3);
        S56SUBSelasLambda = subs(S56elasLambda,d,deltasub);
    end
end 

%% SENSITIVITY OF LAMBDA WITH RESPECT TO PROPORTION MIGRATING RATES
% MIGRATION HAPPENS BEFORE DEMOGRAPHY 
% there are 4 proportion rates for both stage 4 & 5, but within each stage,
% the proportion rates can only be perturbed in pairs, as the columns of
% P^i_k need to sum to 1. Thus there are only two available perturbations
% for each stage, and so 4 perturbations in total. 

% the sensitivity of \bC(\sP) is 0 for all as there is a specified path in
% season 6 (Z_{2,6}=0). season 6 is also when we are performing the
% perturbations. 

% sensitivity formula for perturbations to P^4_6
P46sensLambda = sym(zeros(1,2)); 
% elasticity formula for perturbations to P^4_6 life rates
P46elasLambda = sym(zeros(1,2)); 
% sensitivity formula for perturbations to P^5_6 
P56sensLambda = sym(zeros(1,2)); 
% elasticity formula for perturbations to P^5_6 life rates
P56elasLambda = sym(zeros(1,2)); 

for pp = 1:2
    syms d % act as delta

    if pp == 1 
    % perturbation is to P^4_6(3,2) and P^4_6(4,2) and stage 5 equivalents
        ploc = [0 0 0 0; 0 0 0 0; 0 d 0 0; 0 -d 0 0];
        % 1% of life rate being positively perturbed for elasticity substitutions
        deltasubs = 0.556*0.01;
    else
    % perturbation is to P^4_6(3,3) and P^4_6(4,3) and stage 5 equivalents
        ploc = [0 0 0 0; 0 0 0 0; 0 0 d 0; 0 0 -d 0];
        % 1% of life rate being positively perturbed for elasticity substitutions
        deltasubs = 0.532*0.01;
    end
    
    Ecii = zeros(c); Ecii(4,4)=1; %matrix E_{c,ii} i=4
    vecEcii = reshape(Ecii,c*c,1); 

    pP46 = P(:,3*n+1:4*n,6) + ploc; % perturbed version of P^4_6
    vecPp46 = reshape(pP46,n*n,1);
    
    dvecA6P4 = (Z2k(:,:,6)*(vecmS6.*(X1*kron(diff(vecPp46),vecEcii)))).*Z3k;
    dvecAAP4 = kronYY*dvecA6P4;
    
    P46sensLambda(:,pp) = (kron(W',V')/(V'*W))*dvecAAP4;
    P46elasLambda(:,pp) = d/lambda*P46sensLambda(:,pp);
    P46SUBSelasLambda = subs(P46elasLambda,d,deltasub);

    Ecii = zeros(c); Ecii(5,5)=1; %matrix E_{c,ii} i=5
    vecEcii = reshape(Ecii,c*c,1); 
        
    pP56 = P(:,4*n+1:5*n,6) + ploc; % perturbed version of P^5_6
    vecPp56 = reshape(pP56,n*n,1);
    
    dvecA6P5 = (Z2k(:,:,6)*(vecmS6.*(X1*kron(diff(vecPp56),vecEcii)))).*Z3k;
    dvecAAP5 = kronYY*dvecA6P5;
    
    P56sensLambda(:,pp) = (kron(W',V')/(V'*W))*dvecAAP5;
    P56elasLambda(:,pp) = d/lambda*P56sensLambda(:,pp);
    P56SUBSelasLambda = subs(P56elasLambda,d,deltasub);
end 


%% CALCULATING BY SPECIFYING DELTA SMALL 
d = [-0.000000000001,0.000000000001];% possible values of delta
%% HOW PERTURBATIONS TO DEMOGRAPHY AFFECT ASYMPTOTICS
% annual cycle matrix following perturbations to D_{3,6} life rates
% D36sA = zeros(c*n,c*n,7); % no path specified for \sA
D_36 = D(:,(3-1)*c+1:3*c,6);
D_46 = D(:,(4-1)*c+1:4*c,6);
D36mD = mD(:,:,6); % to store perturbed version of \mD_6
D46mD = mD(:,:,6);
lambdaD36sA = zeros(2,7);
lambdaD46sA = zeros(2,7);
senslambdaD36 = zeros(1,7);
senslambdaD46 = zeros(1,7);
elaslambdaD36 = zeros(1,7);
elaslambdaD46 = zeros(1,7);
for pp = 1:7  %number of non-zero entries in D_36 and D_46
        if pp==1 % perturbation is to D_j6(1,3)
            for dd = 1:2 
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
            for dd = 1:2 
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
            for dd = 1:2
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
            for dd = 1:2 
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
            for dd = 1:2 
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
            for dd = 1:2 
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
            for dd = 1:2 
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
        
%         if pp == 1 || pp == 2 || pp == 3    % reproductive life rates 
%             figure(1)
%             subplot(1,2,1)
%             plot(d,lambdaD36sA(:,pp),'LineWidth',2)
%             hold on
%             plot(d,lambdaD46sA(:,pp),'--','LineWidth',2)
%         else 
%             figure(1)
%             subplot(1,2,2)
%             plot(d,lambdaD36sA(:,pp),'LineWidth',2)
%             hold on
%             plot(d,lambdaD46sA(:,pp),'--','LineWidth',2)
%         end 
        
        % estimate of slope at delta = 0 
        senslambdaD36(pp) = (lambdaD36sA(1,pp)-lambdaD36sA(2,pp))/(d(1)-d(2));
        senslambdaD46(pp) = (lambdaD46sA(1,pp)-lambdaD46sA(2,pp))/(d(1)-d(2));

%         elaslambdaD36(pp) = (d/Cptilde(ll,ii))*P56sensCPtildestart(ii,pp,ll);        
end
% 
% figure(1)
% subplot(1,2,1)
% legend('delta1D36','delta1D46','delta2D36','delta2D46','delta3D36','delta3D46')
% hold on
% subplot(1,2,2)
% legend('delta4D36','delta4D46','delta5D36','delta5D46','delta6D36','delta6D46','delta7D36','delta7D46')

%% HOW PERTRUBATIONS TO MOVEMENT SURVIVAL AFFECT ASYMPTOTICS 
% annual cycle matrix following perturbations to S^i_6, where i \in {4,5}
S_46 = S(:,(4-1)*n+1:4*n,6);
S_56 = S(:,(5-1)*n+1:5*n,6);
S46mS = mS(:,:,6); % to store perturbed version of \mS_6
S56mS = mS(:,:,6);
lambdaS46sA = zeros(2,3); % there are 3 possible entries to perturb
lambdaS56sA = zeros(2,3);
senslambdaS46 = zeros(1,3);
senslambdaS56 = zeros(1,3);

for pp = 1:3  %number of entries in S^4_6 and S^5_6 subject to perturbation
        if pp==1 % perturbation is to S^i_6(3,2)
            for dd = 1:2
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
            for dd = 1:2
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
            for dd = 1:2
                pS46mS = S46mS;
                % perturb the (i,i)-th entry of (4,3)_th block of mS_6. i=4
                pS46mS((4-1)*c+4,(3-1)*c+4) = S46mS((4-1)*c+4,(3-1)*c+4) + d(dd); 
                S46sA = Y127*mD(:,:,6)*(pS46mS.*mP(:,:,6))*Y51;
                lambdaS46sA(dd,pp) = max(abs(eig(S46sA)));
                
                % perturb the (i,i)-th entry of (4,3)_th block of mS_6. i=5
                pS56mS = S56mS;
                pS56mS((3-1)*c+5,(2-1)*c+5) = S56mS((4-1)*c+5,(3-1)*c+5) + d(dd); 
                S56sA = Y127*mD(:,:,6)*(pS56mS.*mP(:,:,6))*Y51;
                lambdaS56sA(dd,pp) = max(abs(eig(S56sA)));
            end
        end
        
%         figure(2)
%         plot(d,lambdaS46sA(:,pp),'LineWidth',2)
%         hold on
%         plot(d,lambdaS56sA(:,pp),'--','LineWidth',2)
%         hold on
%         
       
        % estimate of slope at delta = 0 
        senslambdaS46(pp) = (lambdaS46sA(1,pp)-lambdaS46sA(2,pp))/(d(1)-d(2));
        senslambdaS56(pp) = (lambdaS56sA(1,pp)-lambdaS56sA(2,pp))/(d(1)-d(2));
        
end

% figure(2)
% legend('delta1S46','delta1S56','delta2S46','delta2S56','delta3S46','delta3S56')

%% HOW PERTRUBATIONS TO MOVEMENT PROPORTION AFFECT ASYMPTOTICS 
% annual cycle matrix following perturbations to P^i_6, where i \in {4,5}
P_46 = P(:,(4-1)*n+1:4*n,6);
P_56 = P(:,(5-1)*n+1:5*n,6);
P46mP = mP(:,:,6); % to store perturbed version of \mP_6
P56mP = mP(:,:,6);
lambdaP46sA = zeros(2,2); %there are 4 possible entries to perturb but
lambdaP56sA = zeros(2,2);  %two rates must be perturbed at once to ensure
senslambdaP46 = zeros(1,2);  %proportions sum to 1, so have 2 deltas
senslambdaP56 = zeros(1,2);

for pp = 1:2  %number of deltas that perturb P^4_6 and P^5_6
        if pp==1 % perturbation is to P^i_6(3,2) and P^i_^(4,2)
            for dd = 1:2
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
            for dd = 1:2
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
                pP56mP((3-1)*c+4,(3-1)*c+5) = P56mP((3-1)*c+5,(3-1)*c+5) + d(dd);
                % -ve perturbation to the (i,i)-th entry of (4,3)_th block
                % of mP_6. i=5
                pP56mP((3-1)*c+4,(3-1)*c+5) = P56mP((3-1)*c+5,(3-1)*c+5) - d(dd);
                P56sA = Y127*mD(:,:,6)*(mS(:,:,6).*pP56mP)*Y51;
                lambdaP56sA(dd,pp) = max(abs(eig(P56sA)));
            end
        end
        
%         figure(3)
%         plot(d,lambdaP46sA(:,pp),'LineWidth',2)
%         hold on
%         plot(d,lambdaP56sA(:,pp),'--','LineWidth',2)
%         hold on
%         
        % estimate of slope at delta = 0 
        senslambdaP46(pp) = (lambdaP46sA(1,pp)-lambdaP46sA(2,pp))/(d(1)-d(2));
        senslambdaP56(pp) = (lambdaP56sA(1,pp)-lambdaP56sA(1,pp))/(d(1)-d(2));
        
end

% figure(3)
% legend('delta1P46','delta1P56','deltaPS46','deltaPS56')


