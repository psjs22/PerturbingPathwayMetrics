function  [A,Ahat,n,c,s,mP,D,P,S,M,mD,mS,mM] = MigModel(mm)
% Compute the seasonal matrices for the migratory models in paper
% mm defines the model number: 
% - mm = 1 -> Partial migrants w/ residents in both habitats 
% - mm = 2 -> Partial migrants w/ residents in overwintering habitat
% - mm = 3 -> Partial migrants w/ residents in breeding habitat
% - mm = 4 -> Full migration 
% - mm = 5 -> Four season hypothetical model 
% - mm = 6 -> monarch butterfly model
% - mm = 7 -> Four season hypothetical model w/ different migration rates
% for different stage classes 
% - mm = 8 -> hypothetical migratory population Wiederholt
% - mm = 9 -> Sample's monarch butterfly model

%% variables that are the same for all models
if mm <= 4  % n,c and s are same for mm in {1,2,3,4}
    n = 2;  % number of patches. counted by ii
    c = 2;  % number of stage classes in population. counted by jj
    s = 2;  % number of seasons. counted by kk
elseif mm == 5 || mm == 7  % four season models
    n = 2;
    c = 2;
    s = 4;
elseif mm == 6 % monarch butterfly model
    n = 4;
    c = 5; 
    s = 12;
elseif mm == 8 % Wiederholt model
    n = 4;
    c = 2; 
    s = 2; 
elseif mm == 9 % Sample's monarch butterfly model
    n = 4; 
    c = 1; 
    s = 7; 
end 

%% construct the demographic (D) and movement matrices (P,S & M) 
% demographic matrices
D = zeros(c,c*n,s);     %store s \mD matrices 
% movement matrices
P = zeros(n,c*n,s);     % store s \mP matrices
S = P;     % store s \mS matrices
M = P;     % store s \mM matrices

% constructing the other migratory models
if mm == 1 % Partial migrants w/ residents in both habitats 
    %demographic matrices
    D(:,1:c,1) = [0 0.6665; 0.8 0.9];       %D_{1,1) 
    D(:,1:c,2) = [0.8 0; 0 0.9];
    D(:,1+c:2*c,1) = [0 0.5813; 0.6 0.7];   %D_{2,1}
    D(:,c+1:2*c,2) = [0.6 0; 0 0.7];
    % movement matrices
    P(:,1:n,1) = [0.6 0.4; 0.4 0.6];    
    P(:,n+1:2*n,1) = P(:,1:n,1);
    P(:,:,2) = P(:,:,1); 
    %
    S(:,1:n,1) = [1 0.8; 0.8 1];    
    S(:,n+1:2*n,1) = S(:,1:n,1);
    S(:,:,2) = S(:,:,1); 
    %
    M(:,1:n,1) = P(:,1:n,1).*S(:,1:n,1); %[0.6 0.32; 0.32 0.6];    
    M(:,n+1:2*n,1) = M(:,1:n,1);
    M(:,:,2) = M(:,:,1); 
elseif mm == 2 % Partial migrants w/ residents in overwintering habitat
    % demographic matrices
    D(:,1:c,1) = [0 0.6665; 0.8 0.9];       %D_{1,1) 
    D(:,1:c,2) = eye(c); 
    D(:,1+c:2*c,1) = [0 0.5813; 0.6 0.7];   %D_{2,1}
    D(:,c+1:2*c,2) = [0.6 0; 0 0.7];
    % movement matrices
    M(:,1:n,1) = [0 0.32; 0 0.6]; 
    M(:,n+1:2*n,1) = M(:,1:n,1);
    M(:,1:n,2) = [0 0; 0.8 1]; 
    M(:,n+1:2*n,2) = M(:,1:n,2);
elseif mm ==3 % Partial migrants w/ residents in breeding habitat
    % demographic matrices
    D(:,1:c,1) = [0 0.6665; 0.8 0.9];       %D_{1,1) 
    D(:,1:c,2) = [0.8 0; 0 0.9];
    D(:,1+c:2*c,1) = eye(c); 
    D(:,c+1:2*c,2) = [0.6 0; 0 0.7];
    % movement matrices
    M(:,1:n,1) = [1 0; 0.8 0];
    M(:,n+1:2*n,1) = M(:,1:n,1);
    M(:,1:n,2) = [0.6 0.32; 0 0];
    M(:,n+1:2*n,2) = M(:,1:n,2);
elseif mm == 4 % Full migration
    % demographic matrices
    D(:,1:c,1) = [0 0.6665; 0.8 0.9];       %D_{1,1) 
    D(:,1:c,2) = eye(c);
    D(:,1+c:2*c,1) = eye(c); 
    D(:,c+1:2*c,2) = [0.6 0; 0 0.7];
    % movement matrices
    M(:,1:n,1) = [0 0; 0.8 0]; 
    M(:,n+1:2*n,1) = M(:,1:n,1);
    M(:,1:n,2) = [0 0.8; 0 0]; 
    M(:,n+1:2*n,2) = M(:,1:n,2);
elseif mm == 5 % Four season model 
    % demographic matrices 
    D(:,1:c,1) = [0 0.6665; 0.8 0.9];
    D(:,1:c,2) = [0.8 0; 0 0.9];
    D(:,1:c,3) = [0.7 0; 0 0.8];
    D(:,1:c,4) = [0.6 0; 0 0.7];
    D(:,1+c:2*c,1) = [0 0.5813; 0.6 0.7]; 
    D(:,c+1:2*c,2) = [0.6 0; 0 0.7];
    D(:,1+c:2*c,3) = [0.7 0; 0 0.8]; 
    D(:,c+1:2*c,4) = [0.8 0; 0 0.9];
    % movement matrices
    P(:,1:n,1) = [0.6 0.4; 0.4 0.6];    
    P(:,n+1:2*n,1) = P(:,1:n,1);
    P(:,:,2) = P(:,:,1); 
    P(:,:,3) = P(:,:,1); 
    P(:,:,4) = P(:,:,1); 
    %
    S(:,1:n,1) = [1 0.8; 0.8 1];    
    S(:,n+1:2*n,1) = S(:,1:n,1);
    S(:,:,2) = S(:,:,1); 
    S(:,:,3) = S(:,:,1); 
    S(:,:,4) = S(:,:,1); 
    %
    M(:,1:n,1) = P(:,1:n,1).*S(:,1:n,1); %[0.6 0.32; 0.32 0.6];    
    M(:,n+1:2*n,1) = M(:,1:n,1);
    M(:,:,2) = M(:,:,1); 
    M(:,:,3) = M(:,:,1); 
    M(:,:,4) = M(:,:,1); 
elseif mm == 6 % monarch butterfly model
    % demographic matrices
    s1 = [0.041601 0.041601 0.041601 0.041601;... %january
        0.041601 0.041601 0.041601 0.041601; ... %february
        0.041601 0.041601 0.041601 0.041601; ... %march
        0.041601 0.041601 0.041601 0.041601; ... %april
        0.041601 0.041509 0.041601 0.041601; ... %may
        0.041601 0.034856 0.041598 0.041601; ... %june
        0.041601 0.041061 0.041007 0.041336; ... %july
        0.041601 0.041601 0.01615 0.022597; ... %august
        0.041601 0.041601 0.037292 0.035425; ... %september
        0.041601 0.041025 7.64283*10^(-5) 9.49302*10^(-5); ... %october
        0.041601 0.028151859 0.041601 0.041601; ... %november
        0.041601 0.041601 0.041601 0.041601]; %december
    s2 = 0.9896; s3 = 0.9896; s4 = 0.56; s5 = 0.17; % survival rates 
    f4 = 268; f5 = 89; % fecundity rates, average number of eggs per stem
    d = [1 1 1 1; 1 1 1 1; 0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0;...
        1 0.25 0.25 0.25; 1 0.5 0.5 0.5; 1 1 1 1; 1 1 1 1; ...
        1 1 1 1]; % proportion entering diapause each season
    e = [0 0 0 0; 0 0 0 0; 1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1;...
        0 0.5 0.5 0.5; 0 0.25 0.25 0.25; 0 0 0 0; 0 0 0 0; ...
        0 0 0 0]; % proportion becoming reprodcutive each season
    for kk = 1:s
        D(:,1:c,kk) = [0 0 0 0 0; % no breeding
            s1(kk,1)*d(kk,1) 0 0 0 0; 0 s2*(1-e(kk,1)) s3*(1-e(kk,1)) 0 0;
            s1(kk,1)*(1-d(kk,1)) 0 0 0 0; 0 0 s4*e(kk,1) s4 s5]; % Mexico
        D(:,c+1:2*c,kk) = [0 0 f4*(1-d(kk,2))*e(kk,2) f4 f5;
            s1(kk,2)*d(kk,2) 0 0 0 0; 0 s2*(1-e(kk,2)) s3*(1-e(kk,2)) 0 0;
            s1(kk,2)*(1-d(kk,2)) 0 0 0 0; 0 0 s4*e(kk,2) s4 s5]; % South
        D(:,2*c+1:3*c,kk) = [0 0 f4*(1-d(kk,3))*e(kk,3) f4 f5;
            s1(kk,3)*d(kk,3) 0 0 0 0; 0 s2*(1-e(kk,3)) s3*(1-e(kk,3)) 0 0;
            s1(kk,3)*(1-d(kk,3)) 0 0 0 0; 0 0 s4*e(kk,3) s4 s5]; % Central
        D(:,3*c+1:4*c,kk) = [0 0 f4*(1-d(kk,4))*e(kk,4) f4 f5; 
            s1(kk,4)*d(kk,4) 0 0 0 0; 0 s2*(1-e(kk,4)) s3*(1-e(kk,4)) 0 0;
            s1(kk,4)*(1-d(kk,4)) 0 0 0 0; 0 0 s4*e(kk,4) s4 s5]; % North
    end 
    % movement matrices
    % proportion migrating 
    P(:,n+1:2*n,1) = [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; % stage 2
    P(:,n+1:2*n,2) = P(:,n+1:2*n,1); 
    P(:,n+1:2*n,3) = P(:,n+1:2*n,1);
    P(:,n+1:2*n,4) = zeros(n);
    P(:,n+1:2*n,5) = zeros(n);
    P(:,n+1:2*n,6) = zeros(n);
    P(:,n+1:2*n,7) = zeros(n);
    P(:,n+1:2*n,8) = zeros(n); 
    P(:,n+1:2*n,9) = [0 0 0 0; 0 0 0 0.58; 0 0 0 0.42; 0 0 0 0]; 
    P(:,n+1:2*n,10) = [0 0 0.45 0; 0 0 0.55 0; 0 0 0 0; 0 0 0 0]; %guess
    P(:,n+1:2*n,11) = [1 1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
    P(:,n+1:2*n,12) = P(:,n+1:2*n,1); 
    %
    P(:,3*n+1:4*n,1) = zeros(n);  % stage 4
    P(:,3*n+1:4*n,2) = zeros(n); 
    P(:,3*n+1:4*n,3) = zeros(n); 
    P(:,3*n+1:4*n,4) = [0 0 0 0; 1 0 0 0; 0 0 0 0; 0 0 0 0];
    P(:,3*n+1:4*n,5) = [0 0 0 0; 0 0.62 0 0; 0 0.38 0 0; 0 0 0 0]; 
    P(:,3*n+1:4*n,6) = [0 0 0 0; 0 0 0 0; 0 0.556 0.532 0; 0 0.444 0.468 0]; 
    P(:,3*n+1:4*n,7) = [0 0 0 0; 0 0 0 0; 0 0 0.426 0.183; 0 0 0.574 0.817];
    P(:,3*n+1:4*n,8) = [0 0 0 0; 0 0 0 0; 0 0 0.7 0.278; 0 0 0.3 0.722];
    P(:,3*n+1:4*n,9) = [0 0 0 0; 0 0 0.55 0.58; 0 0 0.45 0.42; 0 0 0 0]; 
    P(:,3*n+1:4*n,10) = [0 0 0 0; 0 1 1 0; 0 0 0 0; 0 0 0 0]; 
    P(:,3*n+1:4*n,11) = zeros(n); 
    P(:,3*n+1:4*n,12) = zeros(n);
    % 
    % survival of migration
    S(:,n+1:2*n,1) = [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; % stage 2
    S(:,n+1:2*n,2) = S(:,n+1:2*n,1); 
    S(:,n+1:2*n,3) = S(:,n+1:2*n,1);
    S(:,n+1:2*n,4) = zeros(n);
    S(:,n+1:2*n,5) = zeros(n);
    S(:,n+1:2*n,6) = zeros(n);
    S(:,n+1:2*n,7) = zeros(n);
    S(:,n+1:2*n,8) = zeros(n); 
    S(:,n+1:2*n,9) = [0 0 0 0; 0 0 0 0.5; 0 0 0 0.5; 0 0 0 0]; 
    S(:,n+1:2*n,10) = [0 0 0.196 0; 0 0 0.567 0; 0 0 0 0; 0 0 0 0]; %guess
    S(:,n+1:2*n,11) = [1 0.69 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
    S(:,n+1:2*n,12) = S(:,n+1:2*n,1); 
    %
    S(:,3*n+1:4*n,1) = zeros(n);  % stage 4
    S(:,3*n+1:4*n,2) = zeros(n); 
    S(:,3*n+1:4*n,3) = zeros(n); 
    S(:,3*n+1:4*n,4) = [0 0 0 0; 0.517 0 0 0; 0 0 0 0; 0 0 0 0];
    S(:,3*n+1:4*n,5) = [0 0 0 0; 0 1 0 0; 0 0.733 0 0; 0 0 0 0]; 
    S(:,3*n+1:4*n,6) = [0 0 0 0; 0 0 0 0; 0 0.733 1 0; 0 0.544 0.742 0]; 
    S(:,3*n+1:4*n,7) = [0 0 0 0; 0 0 0 0; 0 0 1 0.5; 0 0 0.742 1];
    S(:,3*n+1:4*n,8) = S(:,3*n+1:4*n,7);
    S(:,3*n+1:4*n,9) = [0 0 0 0; 0 0 0.567 0.5; 0 0 1 0.5; 0 0 0 0]; 
    S(:,3*n+1:4*n,10) = [0 0 0 0; 0 1 0.567 0; 0 0 0 0; 0 0 0 0]; 
    S(:,3*n+1:4*n,11) = zeros(n); 
    S(:,3*n+1:4*n,12) = zeros(n);
    for kk = 1:s 
    % proportion migrating
        P(:,1:n,kk) = [0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]; % stage 1 migration is the same for all kk
        P(:,2*n+1:3*n,kk) = P(:,n+1:2*n,kk); % stage 2&3 same for all kk 
        P(:,4*n+1:5*n,kk) = P(:,3*n+1:4*n,kk); % stage 4&5 same for all kk 
    % survival of migration 
        S(:,1:n,kk) = P(:,1:n,kk); % stage 1 survival is the same for all kk 
        S(:,2*n+1:3*n,kk) = S(:,n+1:2*n,kk); % stage 2&3 same for all kk
        S(:,4*n+1:5*n,kk) = S(:,3*n+1:4*n,kk); % stage 4&5 same for all kk 
    % movement matrices with proportion and survival of migration    
        M(:,1:n,kk) = P(:,1:n,kk).*S(:,1:n,kk); % stage 1
        M(:,n+1:2*n,kk) = P(:,n+1:2*n,kk).*S(:,n+1:2*n,kk); % stage 2
        M(:,2*n+1:3*n,kk) = P(:,2*n+1:3*n,kk).*S(:,2*n+1:3*n,kk); % stage 3
        M(:,3*n+1:4*n,kk) = P(:,3*n+1:4*n,kk).*S(:,3*n+1:4*n,kk); % stage 4
        M(:,4*n+1:5*n,kk) = P(:,4*n+1:5*n,kk).*S(:,4*n+1:5*n,kk); % stage 5
    end
    elseif mm == 7 % Four season model - different migration for stages  
    % demographic matrices 
    D(:,1:c,1) = [0 0.6665; 0.8 0.9];
    D(:,1:c,2) = [0.8 0; 0 0.9];
    D(:,1:c,3) = [0.7 0; 0 0.8];
    D(:,1:c,4) = [0.6 0; 0 0.7];
    D(:,1+c:2*c,1) = [0 0.5813; 0.6 0.7]; 
    D(:,c+1:2*c,2) = [0.6 0; 0 0.7];
    D(:,1+c:2*c,3) = [0.7 0; 0 0.8]; 
    D(:,c+1:2*c,4) = [0.8 0; 0 0.9];
    % movement matrices
    P(:,1:n,1) = [0.6 0.4; 0.4 0.6];    
    P(:,n+1:2*n,1) = [0.3 0.7; 0.3 0.7];
    P(:,:,2) = P(:,:,1); 
    P(:,:,3) = P(:,:,1); 
    P(:,:,4) = P(:,:,1); 
    %
    S(:,1:n,1) = [1 0.8; 0.8 1];    
    S(:,n+1:2*n,1) = S(:,1:n,1);
    S(:,:,2) = S(:,:,1); 
    S(:,:,3) = S(:,:,1); 
    S(:,:,4) = S(:,:,1); 
    %
    M(:,1:n,1) = P(:,1:n,1).*S(:,1:n,1); %[0.6 0.32; 0.32 0.6];    
    M(:,n+1:2*n,1) = M(:,1:n,1);
    M(:,:,2) = M(:,:,1); 
    M(:,:,3) = M(:,:,1); 
    M(:,:,4) = M(:,:,1); 
    %
elseif mm == 8 % hypothetical migratory population Wiederholt
    %demographic matrices
    D(:,1:c,1) = [0 0.96*0.7; 0.96 0.96];       %D_{1,1) 
    D(:,1+c:2*c,1) = [0 0.96*0.55; 0.96 0.96];   %D_{2,1}
    D(:,2*c+1:3*c,2) = [0.98 0; 0 0.98]; %D_{3,2}
    D(:,3*c+1:4*c,2) = [0.86 0; 0 0.86]; %D_{4,2}
    % movement matrices
    P(:,1:n,1) = [0 0 0 0; 0 0 0 0; 0.5 0.5 0 0; 0.5 0.5 0 0];    
    P(:,n+1:2*n,1) = P(:,1:n,1); % stage 2 same as stage 1
    P(:,1:n,2) = [0 0 0.5 0.5; 0 0 0.5 0.5; 0 0 0 0; 0 0 0 0];
    P(:,n+1:2*n,2) = P(:,1:n,2); % stage 2 same as stage 1
    %
    S(:,1:n,1) = [1 0 0 0; 0 1 0 0; 0.8 0.8 1 0; 0.8 0.8 0 1];
    S(:,n+1:2*n,1) = [1 0 0 0; 0 1 0 0; 0.85 0.85 1 0; 0.85 0.85 0 1];
    S(:,1:n,2) = [1 0 0.85 0.85; 0 1 0.85 0.85; 0 0 1 0; 0 0 0 1]; 
    S(:,n+1:2*n,2) = S(:,1:n,2); % stage 2 same as stage 1
    %
    M(:,1:n,1) = P(:,1:n,1).*S(:,1:n,1); 
    M(:,n+1:2*n,1) = P(:,n+1:2*n,1).*S(:,n+1:2*n,1); 
    M(:,1:n,2) = P(:,1:n,2).*S(:,1:n,2); 
    M(:,n+1:2*n,2) = P(:,n+1:2*n,2).*S(:,n+1:2*n,2); 
elseif mm == 9 % Sample's monarch butterfly model
    %demographic matrices
    D(:,1:c,1) = 0.9392001; % Mexico  winter
    D(:,1+c:2*c,1) = 0; % South
    D(:,1+2*c:3*c,1) = 0; % Central 
    D(:,1+3*c:4*c,1) = 0; % North 
    %
    for kk = 2:s
        D(:,1:c,kk) = 0; % Mexico  april
        D(:,1+c:2*c,kk) = 0.3083652; % South
        D(:,1+2*c:3*c,kk) = 0.3083652; % Central 
        D(:,1+3*c:4*c,kk) = 0.3083652; % North 
    end 
    % movement matrices 
    % proportion migrating
    P(:,1:n,1) = [0 1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; %winter 
    P(:,1:n,2) = [0 0 0 0; 0 0.62 0.38 0; 0 0 0 0; 0 0 0 0]; % april 
    P(:,1:n,3) = [0 0 0 0; 0 0 0.556 0.444; 0 0 0.532 0.468; 0 0 0 0]; %may
    P(:,1:n,4) = [0 0 0 0; 0 0 0 0; 0 0 0.426 0.574; 0 0 0.183 0.817]; %jun
    P(:,1:n,5) = [0 0 0 0; 0 0 0 0; 0 0 0.7 0.3; 0 0 0.278 0.722]; %july
    P(:,1:n,6) = [0 0 0 0; 0 0 0 0; 0 0.55 0.45 0; 0 0.58 0.42 0]; %augus
    P(:,1:n,7) = [0 0 0 0; 1 0 0 0; 1 0 0 0; 0 0 0 0]; %september
    % survival of migration
    S(:,1:n,1) = [1 0.517 0.196 0; 0.517 1 0.733 0.544; ...
        0.196 0.733 1 0.742; 0 0.544 0.742 1]; %winter 
    S(:,1:n,2) = [1 0.517 0.196 0; 0.517 1 0.733 0.544; ...
        0.196 0.733 1 0.742; 0 0.544 0.742 1]; % april 
    S(:,1:n,3) = [1 0.517 0.196 0; 0.517 1 0.733 0.544; ...
        0.196 0.733 1 0.742; 0 0.544 0.742 1]; %may
    S(:,1:n,4) = [1 0.517 0.196 0; 0.517 1 0.733 0.544; ...
        0.196 0.733 1 0.742; 0 0.544 0.742 1]; %jun
    S(:,1:n,5) = [1 0.517 0.196 0; 0.517 1 0.733 0.544; ...
        0.196 0.733 1 0.742; 0 0.544 0.742 1]; %july
    S(:,1:n,6) = [1 0.517 0.196 0; 0.69 1 0.733 0.544; ...
        0.196 0.733 1 0.742; 0 0.544 0.742 1]; %august
    S(:,1:n,7) = [1 0 0 0; 0.69 1 0 0; 0.567 0 1 0; 0.5 0 0 1]; %september
    %
    for kk = 1:s
        M(:,1:n,kk) = P(:,1:n,kk).*S(:,1:n,kk); 
    end
end

%% construct seasonal matrices, A_k
mD = zeros(c*n,c*n,s);
for kk = 1:s
    for jj = 1:n
    E = zeros(n); E(jj,jj) = 1;     
    mD(:,:,kk) = mD(:,:,kk) + kron(E,D(:,1+c*(jj-1):c*jj,kk));
    end 
end

mP = zeros(c*n,c*n,s);
mS = mP; 
mM = mP;
for kk = 1:s
    for ii = 1:c 
    E = zeros(c); E(ii,ii) = 1;  
    mP(:,:,kk) = mP(:,:,kk) + kron(P(:,1+n*(ii-1):n*ii,kk),E);
    mS(:,:,kk) = mS(:,:,kk) + kron(S(:,1+n*(ii-1):n*ii,kk),E);
    mM(:,:,kk) = mM(:,:,kk) + kron(M(:,1+n*(ii-1):n*ii,kk),E);
    end
end

A = zeros(c*n,c*n,s);     % store s A matrices
Ahat = A;                 % store s Ahat matrices 
for kk = 1:s
    A(:,:,kk) = mD(:,:,kk)*mM(:,:,kk); % movement happens first
%     A(:,:,kk) = mM(:,:,kk)*mD(:,:,kk); % demography happens first

    Ahat(:,:,kk) = mD(:,:,kk)*mS(:,:,kk); % movement happens first 
%     Ahat(:,:,kk) = mD(:,:,kk)*mS(:,:,kk): % demography happens first 
end
    
        
