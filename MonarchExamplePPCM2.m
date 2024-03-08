%%%%%%%%%%%%%%%%% PERTURBING PATHWAY CONTRIBUTION METRICS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% MONARCH BUTTERFLY EXAMPLE 2 %%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
% Compute matrices for unperturbed model mm in {1,2,3,4,5,6,7,8,9} 
mm = 6; % specify migratory model
[A,Ahat,n,c,s,mP,D,P,S,M,mD,mS,mM] = MigModel(mm);

% SPECIFYING A PATH IN JUST ONE SEASON --- RUN THROUGH ALL OPTIONS
% unperturbed contribution metrics given in Smith et al. 2022
Phi = 0; 
PATH = zeros(n^(length(Phi)+1),s+1,s);
CP = zeros(size(PATH,1),c*n,s);
CPtilde = zeros(size(PATH,1),c*n,s);
CpAvOS = zeros(n,size(PATH,1),s); 
CpAvOH = zeros(c,size(PATH,1),s); 

% COMPUTING UNPERTURBED CONTRIBUTION METRICS
for kk = 1:s
    Phi = kk; % only season kk has a specified path 
    PATH(:,:,kk) = DistinctPaths(n,s,Phi);  % compute all distinct paths
    
    %compute both types of pathway metric and averages
    [CP(:,:,kk),CPtilde(:,:,kk)] = ...
        PathwayMetrics(c,n,s,A,Ahat,Phi,PATH(:,:,kk),mP); 
    
    % compute pathway metrics averaged over stage from Cptilde
    CpAvOS(:,:,kk) = AvOverStage(c,n,x,PATH(:,:,kk),CPtilde(:,:,kk)); 
    
    % compute pathway metrics averaged over habitat from Cptilde
    CpAvOH(:,:,kk) = AvOverHabitat(c,n,x,PATH(:,:,kk),CPtilde(:,:,kk));
    
    % when kk = 1 can compute habitat metrics from CP
    if kk == 1
        ChP = Cp2Ch(c,n,s,mP,Phi,PATH(:,:,kk),CP(:,:,kk));
    end 
end

%% Applying perturbations of threats
% to store perturbed model matrices
pM = M;     pmM = zeros(c*n,c*n,s);
pS = S;     pmS = zeros(c*n,c*n,s);
pD = D;     pmD = zeros(c*n,c*n,s);
pP = P;     pmP = zeros(c*n,c*n,s);
pAhat = zeros(c*n,c*n,s);
pA = zeros(c*n,c*n,s);

pCP = zeros(size(PATH,1),c*n,s);
pCPtilde = zeros(size(PATH,1),c*n,s);
pCpAvOS = zeros(n,size(PATH,1),s); 
pCpAvOH = zeros(c,size(PATH,1),s); 

for kk = 1:s
    %%% THREATS
    
    % High heat in Southern US limits ability to reproduce successfully 
    % during summer months (Malcolom et al. 1987; Pocius et al. 2021)
    % Southern habitat is habitat 2. Summer months June, July and August =>  
    % D_{2,k} fecundity perturbed for k \in {6,7,8} 
    if kk == 6 || kk == 7 || kk == 8 
        pD(1,6:10,kk) = D(1,6:10,kk)*0.99; % 1% decrease in fecundity
    end
    
    % Overwintering ground (habitat 1) reduced due to illegal logging 
    % (Brower 2012;Vidal 2014). Reduce survival in D_{1,k} 
    pD(1:5,1:5,kk) = D(1:5,1:5,kk)*0.95; % 5% deccrease in survival
    
    % Reduced availability of milkweed (Brower 2012; Pleasants 2013; 
    % Flockhart 2015). Survival of larvae declines with decline of milkweed
    % in habitats 2, 3 and 4
    pD(2,6:20,kk) = D(2,6:20,kk)*0.98;    
    
    %%% UPDATE MATRICES FOLLOWING PERTURBATIONS 
    % update demography in mD matrices
    for jj = 1:n
        E = zeros(n); E(jj,jj) = 1;     
        pmD(:,:,kk) = pmD(:,:,kk) + kron(E,pD(:,1+c*(jj-1):c*jj,kk));
    end 
    
    % update movement matrces mP, mS and mM
    for ii = 1:c 
        E = zeros(c); E(ii,ii) = 1;  
        pmP(:,:,kk) = pmP(:,:,kk) + kron(pP(:,1+n*(ii-1):n*ii,kk),E);
        pmS(:,:,kk) = pmS(:,:,kk) + kron(pS(:,1+n*(ii-1):n*ii,kk),E);
        pmM(:,:,kk) = pmM(:,:,kk) + kron(pM(:,1+n*(ii-1):n*ii,kk),E);
    end
    
    % update A and Ahat matrices
    pA(:,:,kk) = pmD(:,:,kk)*pmM(:,:,kk); % movement happens first
    pAhat(:,:,kk) = pmD(:,:,kk)*pmS(:,:,kk); % movement happens first 
end

% compute perturbed contribution metrics perturbed with threats
for  kk = 1:s
    % Run through all options of specifying path in summer months. 
    Phi = kk; % only season kk has a specified path 
    PATH(:,:,kk) = DistinctPaths(n,s,Phi);  % compute all distinct paths
    
    %compute both types of pathway metric and averages
    [pCP(:,:,kk),pCPtilde(:,:,kk)] = ...
        PathwayMetrics(c,n,s,pA,pAhat,Phi,PATH(:,:,kk),pmP); 
    
    % compute pathway metrics averaged over stage from Cptilde
    pCpAvOS(:,:,kk) = AvOverStage(c,n,x,PATH(:,:,kk),pCPtilde(:,:,kk)); 
    
    % compute pathway metrics averaged over habitat from Cptilde
    pCpAvOH(:,:,kk) = AvOverHabitat(c,n,x,PATH(:,:,kk),pCPtilde(:,:,kk));
    
    % when kk = 1 can compute habitat metrics from CP
    if kk == 1
        pChP = Cp2Ch(c,n,s,pmP,Phi,PATH(:,:,kk),pCP(:,:,kk));
    end 
end
% Change in contribution metrics following threats
CP_change = pCP - CP;
CPtilde_change = pCPtilde - CPtilde;
CpAvOS_change = pCpAvOS - CpAvOS;
CpAvOH_change = pCpAvOS - CpAvOS;
ChP_change = pChP - ChP;

%% Applying perturbations of threats AND conservation 
% same threats are applied to the model as above
p2M = pM;     p2mM = zeros(c*n,c*n,s);
p2S = pS;     p2mS = zeros(c*n,c*n,s);
p2D = pD;     p2mD = zeros(c*n,c*n,s);
p2P = pP;     p2mP = zeros(c*n,c*n,s);
p2Ahat = pAhat;
p2A = pA;

p2CP = zeros(size(PATH,1),c*n,s);
p2CPtilde = zeros(size(PATH,1),c*n,s);
p2CpAvOS = zeros(n,size(PATH,1),s); 
p2CpAvOH = zeros(c,size(PATH,1),s); 

for kk = 1:s
    %%% CONSERVATION
    
    % act to combat illegal logging  to improve survival in habitat 1
    pD(1:5,1:5,kk) = D(1:5,1:5,kk)*1.05; % 5% increase in survival
    % still corresponds to an overall decrease of these survival rates
    
    % improve survival of larvae by improving milkweed sources in habitats
    % 2, 3 and 4. 
    p2D(2,6:20,kk) = pD(2,6:20,kk)*1.1; % increase by 10%  
    
    % improve survival along pathways between habitats 2, 3 and 4 due to
    % improved nectar sources and increased milkweed abundance 
    if kk == 5 || kk == 6 || kk == 7 || kk == 8 
        % seasons where stages 4 and 5 use pathways between habitats 2,3&4
        
        for ii = 4:5    % count through stages
            for rr = 1:n    % count through rows of S^i_k
            for cc = 1:n    % count through cols of S^i_k 
                % DO NOT apply perturbation IF: entry is in row 1 of S^i_k, 
                % OR entry is in col 1 of S^i_k (=> cc \in {1,5,11,16}), 
                % OR when entry is on the diagonal of S^i_k 
                % OR when entry (rr,cc) of S^i_k = 1 
                if rr == 1 || cc == 1 || cc == 5 || cc == 11 || cc == 16 ...
                   || rr == cc
                    continue;    % all stages are unperturbed
                else
                 % perturbations to stages 4 and 5, corresponding to
                 % habitats 2, 3 and 4. Increase by 10%
                    p2S(rr,(ii-1)*n+cc,kk) = pS(rr,(ii-1)*n+cc,kk)*1.1;
                end
            end
            end
        end

    elseif kk == 9 || kk == 10
        % seasons where stages 2-5 use pathways between habitats 2,3&4
        
        for ii = 2:5   % count through stages
            for rr = 1:n    % count through rows of S^i_k
            for cc = 1:n    % count through cols of S^i_k 
                % DO NOT apply perturbation IF: entry is in row 1 of S^i_k, 
                % OR entry is in col 1 of S^i_k (=> cc \in {1,5,11,16}), 
                % OR when entry is on the diagonal of S^i_k 
                % OR when entry (rr,cc) of S^i_k = 1 
                if rr == 1 || cc == 1 || cc == 5 || cc == 11 || cc == 16 ...
                   || rr == cc
                    continue;    % all stages are unperturbed
                else
                 % perturbations to stages 2, 3, 4 and 5, corresponding to
                 % habitats 2, 3 and 4 . Increase by 10%
                    p2S(rr,(ii-1)*n+cc,kk) = pS(rr,(ii-1)*n+cc,kk)*1.1;
                end
            end
            end
        end
% 
%     else % No stages are perturbed 
%         continue; % all stages are unperturbed
    end
    
    %%% UPDATE MATRICES FOLLOWING PERTURBATIONS 
    % update M matrices following perturbations to P and S matrices
    p2M = p2S.*p2P;
    % update demography in mD matrices
    for jj = 1:n
        E = zeros(n); E(jj,jj) = 1;     
        p2mD(:,:,kk) = p2mD(:,:,kk) + kron(E,p2D(:,1+c*(jj-1):c*jj,kk));
    end 
    % update movement matrces mP, mS and mM
    for ii = 1:c 
        E = zeros(c); E(ii,ii) = 1;  
        p2mP(:,:,kk) = p2mP(:,:,kk) + kron(p2P(:,1+n*(ii-1):n*ii,kk),E);
        p2mS(:,:,kk) = p2mS(:,:,kk) + kron(p2S(:,1+n*(ii-1):n*ii,kk),E);
        p2mM(:,:,kk) = p2mM(:,:,kk) + kron(p2M(:,1+n*(ii-1):n*ii,kk),E);
    end
    % update A and Ahat matrices
    p2A(:,:,kk) = p2mD(:,:,kk)*p2mM(:,:,kk); % movement happens first
    p2Ahat(:,:,kk) = p2mD(:,:,kk)*p2mS(:,:,kk); % movement happens first   
end

% compute contribution metrics perturbed with threats and conservation
% efforts
for kk = 1:s
    % Run through all options of specifying path in summer months. 
    Phi = kk; % only season kk has a specified path 
    PATH(:,:,kk) = DistinctPaths(n,s,Phi);  % compute all distinct paths
    
    %compute both types of pathway metric and averages
    [p2CP(:,:,kk),p2CPtilde(:,:,kk)] = ...
        PathwayMetrics(c,n,s,p2A,p2Ahat,Phi,PATH(:,:,kk),p2mP); 
    
    % compute pathway metrics averaged over stage from Cptilde
    p2CpAvOS(:,:,kk) = AvOverStage(c,n,x,PATH(:,:,kk),p2CPtilde(:,:,kk)); 
    
    % compute pathway metrics averaged over habitat from Cptilde
    p2CpAvOH(:,:,kk) = AvOverHabitat(c,n,x,PATH(:,:,kk),p2CPtilde(:,:,kk));
    
    % when kk = 1 can compute habitat metrics from CP
    if kk == 1
        p2ChP = Cp2Ch(c,n,s,p2mP,Phi,PATH(:,:,kk),p2CP(:,:,kk));
    end  
end
% Change in contribution metrics following threats and conservation
CP_change2 = p2CP - CP;
CPtilde_change2 = p2CPtilde - CPtilde;
CpAvOS_change2 = p2CpAvOS - CpAvOS;
CpAvOH_change2 = p2CpAvOS - CpAvOS;
ChP_change2 = p2ChP - ChP;

%% Plot C-metrics all on one plot

% condense the C-metrics so that only non-zero rows are stored
% non-zero entries are equal to each other 
cCP = [];
cpCP = [];
cp2CP = [];
cCP_change = [];
cCP_change2 = [];
%
cCPtilde = [];
cpCPtilde = [];
cp2CPtilde = [];
cCPtilde_change = [];
cCPtilde_change2 = [];
%
cCpAvOH = [];
cpCpAvOH = [];
cp2CpAvOH = [];
cCpAvOS = [];
cpCpAvOS = [];
cp2CpAvOS = [];

for kk = 1:s
    for ii = 1:height(CP)
        if any(CP(ii,:,kk)) % there is a non-zero entry in row ii 
            % only cols 2 and 3 are populated and they are equal to each
            % other, so just store col 2
            cCP  = [cCP; CP(ii,2,kk)]; 
            cpCP = [cpCP; pCP(ii,2,kk)];
            cp2CP = [cp2CP; p2CP(ii,2,kk)];
            cCP_change = [cCP_change; CP_change(ii,2,kk)];
            cCP_change2 = [cCP_change2; CP_change2(ii,2,kk)];
            %
            cCPtilde = [cCPtilde; CPtilde(ii,2,kk)];
            cpCPtilde = [cpCPtilde; pCPtilde(ii,2,kk)];
            cp2CPtilde = [cp2CPtilde; p2CPtilde(ii,2,kk)];
            cCPtilde_change = [cCPtilde_change; CPtilde_change(ii,2,kk)];
            cCPtilde_change2 = [cCPtilde_change2; CPtilde_change2(ii,2,kk)];
        end
        
    end
    %
    for ii = 1:height(CpAvOH)
        if any(CpAvOH(ii,:,kk)) % there is a non-zero entry in row ii 
            cCpAvOH = [cCpAvOH; CpAvOH(ii,:,kk)];
            cpCpAvOH = [cpCpAvOH; pCpAvOH(ii,:,kk)];
            cp2CpAvOH = [cp2CpAvOH; p2CpAvOH(ii,:,kk)];
            %
            cCpAvOS = [cCpAvOS; CpAvOS(ii,:,kk)];
            cpCpAvOS = [cpCpAvOS; pCpAvOS(ii,:,kk)];
            cp2CpAvOS = [cp2CpAvOS; p2CpAvOS(ii,:,kk)];
        end
    end                
end

% delete every other row of average Cmetrics matrices
cCpAvOH(2:2:end,:) = [];
cpCpAvOH(2:2:end,:) = [];
cp2CpAvOH(2:2:end,:) = [];
%
cCpAvOS(2:2:end,:) = [];
cpCpAvOS(2:2:end,:) = [];
cp2CpAvOS(2:2:end,:) = [];

xCP = 1:height(cCP);
xCP = xCP+0.33;
xCPtilde = 1:height(cCP);
xCPtilde = xCPtilde +0.67;


figure(1)
plot(0:29,zeros(30),'k')     % x=0
hold on;
% plot the change in subpopulation contribution metrics following threats for each individual pathway
plot(xCP,cCP_change,'v','MarkerEdgeColor',[0 0.4470 0.7410],'MarkerFaceColor',[0 0.4470 0.7410],'MarkerSize',5);
hold on;
% plot the change in subpopulation contribution metrics following threats and conservation for each individual pathway
plot(xCP,cCP_change2,'^','MarkerEdgeColor',[0 0.4470 0.7410],'MarkerFaceColor',[0 0.4470 0.7410],'MarkerSize',5);
hold on;
% plot the change in metapopulation contribution metrics following threats for each individual pathway
plot(xCPtilde,cCPtilde_change,'v','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',5);
hold on;
% plot the change in metapopulation contribution metrics following threats and conservation for each individual pathway
plot(xCPtilde,cCPtilde_change2,'^','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',5);
hold on;
% plot the change in habitat contribution metrics following threats
plot(0.5,ChP_change(2),'v','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5)
hold on;
% plot the change in habitat contribution metrics following threats and conservation
plot(0.5,ChP_change2(2),'^','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5)
hold on;
%
xlim([0 29])
grid on
% grid minor
xticks(0:28)
xticklabels([]);
xlh = xlabel('Pathway');
xlh.Position(2)= xlh.Position(2) -0.12;
ylabel('Change in Contribution Metric')

%% Improvement from conservation actions
cCPcons = cp2CP - cpCP;
cCPtildecons = cp2CPtilde - cpCPtilde; 

%% Reorganise storage of C metrics for plotting in R
% store data related to change in C(P) following threats
CPthreats4R = cell(28,4);
for path = 1:28
    CPthreats4R{path,1} = path; 
    CPthreats4R{path,2} = cCP_change(path); % change in C(P)
    CPthreats4R{path,3} = char(84); % perturbation
    CPthreats4R{path,4} = char(83); % metric
end 
% store data related to change in C(P) following threats and conservation
CPcons4R = cell(28,4);
for path = 1:28
    CPcons4R{path,1} = path;
    CPcons4R{path,2} = cCP_change2(path); % change in C(P)
    CPcons4R{path,3} = char(67); % pert
    CPcons4R{path,4} = char(83); % metric
end 
% store data related to change in Ctilde(P) following threats
CPtildethreats4R = cell(28,4); 
for path = 1:28
    CPtildethreats4R{path,1} = path;
    CPtildethreats4R{path,2} = cCPtilde_change(path); % change in Ctilde(P)
    CPtildethreats4R{path,3} = char(84); % pert
    CPtildethreats4R{path,4} = char(77); % metric
end 
% store data related to change in Ctilde(P) following threats and conservation
CPtildecons4R = cell(28,4); % to 
for path = 1:28
    CPtildecons4R{path,1} = path;
    CPtildecons4R{path,2} = cCPtilde_change2(path); % change in Ctilde(P)
    CPtildecons4R{path,3} = char(67); % pert
    CPtildecons4R{path,4} = char(77); % metric
end 
% store data related to change in habitat Cmetric following threats
Chabitatthreats4R = cell(1,4); 
Chabitatthreats4R{1,1} = 0; 
Chabitatthreats4R{1,2} = ChP_change(2); % value
Chabitatthreats4R{1,3} = char(84); % pert
Chabitatthreats4R{1,4} = char(72); % metric 
% store data related to change in habitat Cmetric following threats and conservation
Chabitatcons4R = cell(1,4);
Chabitatcons4R{1,1} = 0; 
Chabitatcons4R{1,2} = ChP_change2(2); % value
Chabitatcons4R{1,3} = char(67); % pert
Chabitatcons4R{1,4} = char(72); % metric 

% Combine data for all Cmetrics above into one dataframe
Monarch2data4R = cell(28*4+2,5);
for ii = 1:length(Monarch2data4R)
    if ii == 1  % Store Habitat metrics pert by threats
        Monarch2data4R{ii,1} = 0; % pathway 0
        Monarch2data4R{ii,2} = ChP_change(2); % value
        Monarch2data4R{ii,3} = char(84); % pert 
        Monarch2data4R{ii,4} = char(72); % metric
        Monarch2data4R{ii,5} = 0.5; % x value in plot
    elseif ii == 2 % Store Habitat metrics pert by threats and cons
        Monarch2data4R{ii,1} = 0; % pathway 0
        Monarch2data4R{ii,2} = ChP_change2(2); % value
        Monarch2data4R{ii,3} = char(67); % pert 
        Monarch2data4R{ii,4} = char(72); % metric
        Monarch2data4R{ii,5} = 0.5; % x value in plot
    elseif ii >= 3 && ii <= 30 % Store Subpop metrics (CP) pert by threats
        Monarch2data4R{ii,1} = ii-2; % path
        Monarch2data4R{ii,2} = cCP_change(ii-2); % value
        Monarch2data4R{ii,3} = char(84); % pert 
        Monarch2data4R{ii,4} = char(83); % metric
        Monarch2data4R{ii,5} = ii-2+0.3; % x value in plot
    elseif ii >= 31 && ii <= 58 % Store Subpop metrics (CP) pert by threats and cons
        Monarch2data4R{ii,1} = ii-30; % path
        Monarch2data4R{ii,2} = cCP_change2(ii-30); % value
        Monarch2data4R{ii,3} = char(67); % pert 
        Monarch2data4R{ii,4} = char(83); % metric
        Monarch2data4R{ii,5} = ii-30+0.3; % x value in plot
    elseif ii >= 59 && ii <= 86 % Store Metapop metrics (CPtilde) pert by threats
        Monarch2data4R{ii,1} = ii-58; % path
        Monarch2data4R{ii,2} = cCPtilde_change(ii-58); % value
        Monarch2data4R{ii,3} = char(84); % pert 
        Monarch2data4R{ii,4} = char(77); % metric
        Monarch2data4R{ii,5} = ii-58+0.7; % x value in plot
    elseif ii >= 87 && ii <= 114 % Store Metapop metrics (CPtilde) pert by threats and cons
        Monarch2data4R{ii,1} = ii-86; % path
        Monarch2data4R{ii,2} = cCPtilde_change2(ii-86); % value
        Monarch2data4R{ii,3} = char(67); % pert 
        Monarch2data4R{ii,4} = char(77); % metric
        Monarch2data4R{ii,5} = ii-86+0.7; % x value in plot
    end
end

% save data as csv file
writecell(Monarch2data4R,'...\Monarch2data.csv'); % save to same location as code 
