%%%% SENSITIVITY OF CONTRIBUTION METRICS. HYPOTHETICAL EXAMPLE 3.2 %%%%

clear all; close all; 

%%% COMPUTE MATRICES OF THE HYPOTHETICAL MODEL
mm = 1; % specify migratory model
[A,Ahat,n,c,s,mP,D,P,S,M,mD,mS,mM] = MigModel(mm);   

%%% COMPUTE UNPERTURBED PCMs
Phi = [1,2];    % full migratory routes
% compute all distict paths
PATH = DistinctPaths(n,s,Phi);  
%compute pathway metrics
[Cp,Ctilde] = PathwayMetrics(c,n,s,A,Ahat,Phi,PATH,mP);

%% structural matrices and constants that are independent of path
% only a function of \bp in habitat 2 (j=2) and season 1 (k=1)

Z11 = kron(mS(:,:,1)',eye(c*n));   %constant Z_{1,1}
%Z41 = eye(c^2*n^2); %constant Z_{4,1}

K = zeros(c*n); % vec permutation matrix K_{2,2}
for i=1:2
    for j=1:2
        E = zeros(c,n);
        E(i,j) = 1; 
        K = K + kron(E,E');
    end
end

X1 =  kron(kron(eye(n),K),eye(c)); % constant X_1

X4 = kron(eye(c*n),ones(1,c*n)); % constant 

vecE222 = [0;0;0;1];
%vecI2 = [1;0;0;1];

%% SENSITIVITY FORMULA
syms p(d) 2 real  % p is a function of delta (denoted d)
p1_1(d) = D(1,3,1); % (D_{2,1})_{1,1} is unchanged 
p1_2(d) = D(1,4,1); % (D_{2,1})_{1,2} is unchanged 
p2_1(d) = D(2,3,1); % (D_{2,1})_{2,1} is unchanged 
p2_2(d) = D(2,4,1) + d; % (D_{2,1})_{2,2} is a linear function of delta

Dp = subs(p);    % new demog matrix for j=2, k=1 as function of delta
vecDp = reshape(Dp,4,1);   %vectorise Dp

%% Pathway 3 = 1->2->1
EkronJP3 = kron([0 0; 1 0],ones(2)); %E_{2,21} \otimes J_2
Z31P3 = EkronJP3(:);      % constant Z_{3,1} for Path 3
A2P3 = Ahat(:,:,2).*kron([0 1; 0 0],ones(c)); %\bA_2 for Path 3
% sensitivity of subpopulation contribution metric using path 3 to D_{2,1}(2,2)
sensCP3 = X4*kron(eye(4),A2P3)*((Z11*X1*kron(vecE222,diff(vecDp))).*Z31P3);

% ONLY D_{2,1} AFFECTED SO DON'T NEED dvecP1/dp^T part
P1P3 = kron(eye(n),mP(3:4,1:2,1));  %\bP_1 for Path 3
P2P3 = kron(eye(n),mP(1:2,3:4,2));  %\bP_2 for Path 3
Upsilon22P3 = A2P3*P2P3; %\Upsilon^2_2 for Path 3
% sensitivity of metapopulation contribution metric using path 3 to D_{2,1}(2,2)
sensCtildeP3 = X4*kron(eye(c*n),Upsilon22P3)...
    *kron(P1P3',eye(c*n))*((Z11*X1*kron(vecE222,diff(vecDp))).*Z31P3);

%% Pathway 4 = 1->2->2
Z31P4 = Z31P3; %Z_{3,1} is the same as for P3 as path used in k=1 the same
A2P4 = Ahat(:,:,2).*kron([0 0; 0 1],ones(2)); %\bA_2 for Path 4
% sensitivity of subpopulation contribution metric using path 4 to D_{2,1}(2,2)
sensCP4 = X4*kron(eye(4),A2P4)*((Z11*X1*kron(vecE222,diff(vecDp))).*Z31P4);

% ONLY D_{2,1} AFFECTED SO DON'T NEED dvecP1/dp^T part
P1P4 = P1P3; %same pathway used in k=1 as for Path 3
P2P4 = kron(eye(n),mP(3:4,3:4,2)); %\bP_2 for Path 4
Upsilon22P4 = A2P4*P2P4; %\Upsilon for Path 4
% sensitivity of metapopulation contribution metric using path 4 to D_{2,1}(2,2)
sensCtildeP4 = X4*kron(eye(c*n),Upsilon22P4)...
    *kron(P1P4',eye(c*n))*((Z11*X1*kron(vecE222,diff(vecDp))).*Z31P4);

%% Pathway 7 = 2->2->1
EkronJP7 = kron([0 0; 0 1],ones(2)); %E_{2,22} \otimes J_2
Z31P7 = EkronJP7(:); %Z_{3,1} for Path 7
A2P7 = A2P3; %same pathway used in k=2 as for P3
% sensitivity of subpopulation contribution metric using path 7 to D_{2,1}(2,2)
sensCP7 = X4*kron(eye(4),A2P7)*((Z11*X1*kron(vecE222,diff(vecDp))).*Z31P7);

% ONLY D_{2,1} AFFECTED SO DON'T NEED dvecP1/dp^T part
P1P7 = kron(eye(n),mP(3:4,3:4,1)); %\bP_1 for Path 7
P2P7 = P2P3; %same pathway used in k=2 as for Path 3
Upsilon22P7 = A2P7*P2P7; %\Upsilon for Path 7
% sensitivity of metapopulation contribution metric using path 7 to D_{2,1}(2,2)
sensCtildeP7 = X4*kron(eye(c*n),Upsilon22P7)...
    *kron(P1P7',eye(c*n))*((Z11*X1*kron(vecE222,diff(vecDp))).*Z31P7);

%% Pathway 8 = 2->2->2
Z31P8 = Z31P7; %Z_{3,1} is the same as for P7 as path used in k=1 the same
A2P8 = A2P4; %same pathway used in k=2 as for P4
% sensitivity of subpopulation contribution metric using path 8 to D_{2,1}(2,2)
sensCP8 = X4*kron(eye(4),A2P8)*((Z11*X1*kron(vecE222,diff(vecDp))).*Z31P8);

% ONLY D_{2,1} AFFECTED SO DON'T NEED dvecP1/dp^T part
P1P8 = P1P7; %same pathway used in k=1 as for Path 7
P2P8 = P2P4; %same pathway used in k=2 as for Path 4
Upsilon22P8 = A2P8*P2P8; %\Upsilon for Path 8
% sensitivity of metapopulation contribution metric using path 8 to D_{2,1}(2,2)
sensCtildeP8 = X4*kron(eye(c*n),Upsilon22P8)...
    *kron(P1P8',eye(c*n))*((Z11*X1*kron(vecE222,diff(vecDp))).*Z31P8);

%% Habitat Contribution Metrics (no specified paths) 
newZ11 = kron((mS(:,:,1).*mP(:,:,1))',eye(c*n)); %Z11 for no specified path
sensC = X4*kron(eye(c*n),A(:,:,2))*((newZ11*X1*kron(vecE222,diff(vecDp)))...
    .*ones(c^2*n^2,1));

