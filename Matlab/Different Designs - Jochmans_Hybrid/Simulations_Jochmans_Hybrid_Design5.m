%%%% Author: Gabriela Szini
%%%% Date: 25/07/2020
%%%% File based on code in
%%%% Prelimiary/Simulations_Jochmans_Hybrid_correctFE_optimize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DESIGN 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [OUTPUT]=Simulations
% parameters
R  = 500; % # of Monte Carlo replications
n  =   25; % # of nodes 
beta_11 = 1;  % parameter value
beta_12 = 2.5;
beta_21 = 0.8;
beta_22 = 1;
beta_23 = 2;

% initialization
CMLE = cell(R, 1); SECMLE = CMLE ; BETAS = CMLE;  FE_SELECTION = CMLE; ZIJ = CMLE; BETAS_REAL = CMLE;
% Replications
parfor r=1:R, [CMLE{r} SECMLE{r} BETAS{r} FE_SELECTION{r} ZIJ{r} BETAS_REAL{r}] = Simulation(n, beta_11, beta_12, beta_21, beta_22, beta_23); end
% prepare screem output
%TCMLE = (CMLE-ones(R,1)*theta)./SECMLE; aCMLE = abs(TCMLE)>=1.96;
%OUTPUT = [CMLE SECMLE TCMLE aCMLE];
OUTPUT = [CMLE SECMLE BETAS FE_SELECTION ZIJ BETAS_REAL]
save('Simulations_Jochmans_Hybrid_Design5_n25.mat','OUTPUT')

function [cmle se_cmle betas FE_selection_real_estimated zij_real_estimated betas_real] = Simulation(n, beta_11, beta_12, beta_21, beta_22, beta_23)
warning('off','all'); % prevent overparametrization warning for MLE

%%%% Data generation % see paper for details on design choices %%%%

% generating FE
a = repmat(1:(n*n),1)';
monte_fe_i_index = sort(a-n*floor(a/n) + 1); 
monte_fe_j_index = (a-n*floor(a/n) + n + 1); 
monte_fe_indices = [monte_fe_i_index,monte_fe_j_index];
monte_fe_indices = sortrows(monte_fe_indices, [1 2]);
reference = monte_fe_indices(:, 2)-n;
condition = monte_fe_indices(:,1)==reference(:); %changed here
monte_fe_indices(condition,:) = [];
monte_fe_indices = [monte_fe_indices,zeros(size(monte_fe_indices,1),2)]; %here as well
% no fixed effects in observation equation
%fixed_effects_observation = normrnd(0, 1, [(n+n), 1]);
fixed_effects_selection = normrnd(0, 1, [(n+n), 1]);
monte_fe_indices(:,3)= fixed_effects_selection(monte_fe_indices(:,1)) - mean(fixed_effects_selection(monte_fe_indices(:,1)));
monte_fe_indices(:,4)= fixed_effects_selection(monte_fe_indices(:,2)) - mean(fixed_effects_selection(monte_fe_indices(:,2)));
%monte_fe_indices(:,5)= fixed_effects_selection(monte_fe_indices(:,1)) - mean(fixed_effects_selection(monte_fe_indices(:,1)));
%monte_fe_indices(:,6)= fixed_effects_selection(monte_fe_indices(:,2)) - mean(fixed_effects_selection(monte_fe_indices(:,2)));

% generating explanatory variables
number_observations = n*(n-1);
monte_x1 = normrnd(0, 1, [number_observations, 1]) + monte_fe_indices(:,3) + monte_fe_indices(:,4);
monte_x2 = rand(number_observations,1)<=0.5;
monte_x3 = rand(number_observations,1)<=0.5;

% generating error terms
mu = [0 0];
sigma = [1 -0.7; -0.7 1];
normal_errors = mvnrnd(mu,sigma,number_observations);
% for observation equation we keep normal errors, for the selection
% equation transform to logistic
errors_observation = normal_errors(:,1);
errors_selection = icdf('logistic',normcdf(normal_errors(:,2)));

%generating dependent variables
monte_y2 = (beta_21*monte_x1 + beta_22*monte_x2 + beta_23*monte_x3 + monte_fe_indices(:,3) + monte_fe_indices(:,4) + errors_selection)>=0;
monte_y1 = (beta_11*monte_x1 + beta_12*monte_x2 + errors_observation).*(monte_y2>0);
monte_y1(monte_y1==0) = NaN;

%put all in a matrix
monte_dgp = [monte_fe_indices(:,1),monte_fe_indices(:,2),monte_y1,monte_x1, monte_x2,monte_y2,monte_fe_indices(:,3),monte_fe_indices(:,4),monte_x3];

%%%% Conditional logit estimation %%%%

%Put back diagonal elements
add =  NaN(n,size(monte_dgp,2));
add(:,1)= repmat(1:n,1)';
add(:,2)= repmat(1:n,1)'+n;
monte_dgp = sortrows(vertcat(monte_dgp,add),[1 2]);

%create ij variable
ij_monte_dgp = strcat(num2str(monte_dgp(:,1)),num2str(zeros(size(monte_dgp,1),1)),num2str(monte_dgp(:,2)));
ij_monte_dgp = mat2cell(ij_monte_dgp, ones(size(ij_monte_dgp,1),1), size(ij_monte_dgp,2));
ij_monte_dgp = cellfun(@str2num, ij_monte_dgp,'UniformOutput',false);

%create matrices
y_matrix_check_indices = reshape(ij_monte_dgp,n,n)';
monte_y2_matrix = reshape(monte_dgp(:,6),n,n)';
monte_x1_matrix = reshape(monte_dgp(:,4),n,n)';
monte_x2_matrix = reshape(monte_dgp(:,5),n,n)';
monte_x3_matrix = reshape(monte_dgp(:,9),n,n)';

% conditional logit estimation
[cmle se_cmle m_star combinations ss] = ConditionalLogit(monte_y2_matrix,monte_x1_matrix,monte_x2_matrix,monte_x3_matrix); 

%%%% Estimating FE in hybrid approach %%%%
[fe_xi_i_zeta_j_est monte_dgp_subset condition_all] = FEHybridApproach(monte_dgp,combinations,ss,n,cmle);
FE_selection_real_estimated = [repmat(1:(n+n),1)', fixed_effects_selection, fe_xi_i_zeta_j_est];

%%%% Estimating the second stage %%%%

[betas zij_real_estimated betas_real] = SecondStage(monte_dgp_subset,condition_all,fe_xi_i_zeta_j_est,cmle,fixed_effects_selection,beta_21,beta_22,beta_23);

function [fe_xi_i_zeta_j_est monte_dgp_subset condition_all] = FEHybridApproach(monte_dgp,combinations,s,n,cmle)
%removing missings in x2
monte_dgp_subset = monte_dgp;
monte_dgp_subset(any(ismissing(monte_dgp_subset(:,7)),2), :) = []; 

%make ij here for dgp_subset
ij_monte_dgp_subset = strcat(num2str(monte_dgp_subset(:,1)),num2str(zeros(size(monte_dgp_subset,1),1)),num2str(monte_dgp_subset(:,2)));
ij_monte_dgp_subset = mat2cell(ij_monte_dgp_subset, ones(size(ij_monte_dgp_subset,1),1), size(ij_monte_dgp_subset,2));
ij_monte_dgp_subset = cell2mat(cellfun(@str2num, ij_monte_dgp_subset,'UniformOutput',false));

%as need combinations used in conditional MLE, we include this variable
combinations_subset_i1 = combinations((s==1),1);
combinations_subset_j1 = combinations((s==1),3)+n;
combinations_subset = unique([combinations_subset_i1, combinations_subset_j1],'rows');

% make ij here for the combinations
ij_combinations_subset = strcat(num2str(combinations_subset(:,1)),num2str(zeros(size(combinations_subset,1),1)),num2str(combinations_subset(:,2)));
ij_combinations_subset = mat2cell(ij_combinations_subset, ones(size(ij_combinations_subset,1),1), size(ij_combinations_subset,2));
ij_combinations_subset = cell2mat(cellfun(@str2num, ij_combinations_subset,'UniformOutput',false));

% check which combinations are present
monte_dgp_subset_subset_index = ismember(ij_monte_dgp_subset,ij_combinations_subset);
monte_dgp_subset_subset = monte_dgp_subset((monte_dgp_subset_subset_index==1),:);

% check which FE we should have estimated at the end
fe_xi_subset= unique(combinations_subset(:,1));
fe_zeta_subset = unique(combinations_subset(:,2));

% checking conditions for fixed effects
% first for fe_xi
fe_delete_xi = [];
fe_delete_zeta = [];
condition_while = 1;
while condition_while == 1,
    Nsum_xi = [];
    Nsum_zeta = [];
    Ncount_xi = histc(monte_dgp_subset_subset(:,1),fe_xi_subset);
    for i=1:size(fe_xi_subset,1),
        Nsum_xi(i) = sum(monte_dgp_subset_subset(monte_dgp_subset_subset(:,1)==fe_xi_subset(i),6));
    end
    Ncondition1_xi = (Nsum_xi'==Ncount_xi);
    Ncondition2_xi = (Nsum_xi'==0);
    fe_delete_new_xi = vertcat(fe_xi_subset(find(Ncondition1_xi == 1)),fe_xi_subset(find(Ncondition2_xi == 1)));
    fe_delete_xi = vertcat(fe_delete_xi,fe_delete_new_xi);
    fe_xi_subset = fe_xi_subset(ismember(fe_xi_subset,fe_delete_xi)==0); 
    % then for fe_zeta
    Ncount_zeta = histc(monte_dgp_subset_subset(:,2),fe_zeta_subset);
    for i=1:size(fe_zeta_subset,1),
        Nsum_zeta(i) = sum(monte_dgp_subset_subset(monte_dgp_subset_subset(:,2)==fe_zeta_subset(i),6));
    end
    Ncondition1_zeta = (Nsum_zeta'==Ncount_zeta);
    Ncondition2_zeta = (Nsum_zeta'==0);
    fe_delete_new_zeta = vertcat(fe_zeta_subset(find(Ncondition1_zeta == 1)),fe_zeta_subset(find(Ncondition2_zeta == 1)));
    fe_delete_zeta = vertcat(fe_delete_zeta,fe_delete_new_zeta);
    fe_zeta_subset = fe_zeta_subset(ismember(fe_zeta_subset,fe_delete_zeta)==0) ;
    conditions = (not(ismember(monte_dgp_subset_subset(:,1),fe_delete_xi))).*(not(ismember(monte_dgp_subset_subset(:,2),fe_delete_zeta)));
    monte_dgp_subset_subset = monte_dgp_subset_subset((conditions==1),:);
    if isempty(fe_delete_new_xi)==1 && isempty(fe_delete_new_zeta)==1,
        condition_while = 0;
    end
end

% another conditions for excluding in likelihood
condition_all = (not(ismember(monte_dgp_subset(:,1),fe_delete_xi))).*(not(ismember(monte_dgp_subset(:,2),fe_delete_zeta))).*monte_dgp_subset_subset_index;

% initial values for log-likelihood
fe_xi_i_zeta_j_est = ones(n+n,1);

% define likelihood
F_logit = @(e) 1./(1+exp(-e));
LikelihoodHybrid = @(fe_xi_i_zeta_j_est) -sum((monte_dgp_subset(:,6).*(log(F_logit(fe_xi_i_zeta_j_est(monte_dgp_subset(:,1))+fe_xi_i_zeta_j_est(monte_dgp_subset(:,2))+cmle(1)*monte_dgp_subset(:,4)+cmle(2)*monte_dgp_subset(:,5)+cmle(3)*monte_dgp_subset(:,9))))+(1-monte_dgp_subset(:,6)).*(log(1-(F_logit(fe_xi_i_zeta_j_est(monte_dgp_subset(:,1))+fe_xi_i_zeta_j_est(monte_dgp_subset(:,2))+cmle(1)*monte_dgp_subset(:,4)+cmle(2)*monte_dgp_subset(:,5)+cmle(3)*monte_dgp_subset(:,9)))))).*condition_all(:)); 

% use fmincon to maximize
options = optimoptions('fmincon','MaxFunctionEvaluations',100000);
result_fe = fmincon(LikelihoodHybrid,fe_xi_i_zeta_j_est,[],[],[],[],[],[],[],options);
    
% for fixed effects : the ones that are not calculated we should put NaN:
no_fe = ismember(repmat(1:(n+n),1)',vertcat(fe_xi_subset,fe_zeta_subset));
fe_xi_i_zeta_j_est = result_fe;
for i=1:size(no_fe),
    if no_fe(i) == 0
        fe_xi_i_zeta_j_est(i) = NaN;
    end
end
   
function [betas zij_real_estimated betas_real] = SecondStage(monte_dgp_subset,condition_all,fe_xi_i_zeta_j_est,cmle,fixed_effects_selection,beta_21,beta_22,beta_23)
% taking only combinations used in Jochmans
monte_dgp_secondstage = monte_dgp_subset((condition_all==1),:);

% create variable zij and getting inverted mills ratio
monte_dgp_secondstage_zij = fe_xi_i_zeta_j_est(monte_dgp_secondstage(:,1)) + fe_xi_i_zeta_j_est(monte_dgp_secondstage(:,2)) + cmle(1)*monte_dgp_secondstage(:,4) + cmle(2)*monte_dgp_secondstage(:,5) + cmle(3)*monte_dgp_secondstage(:,9);
monte_dgp_secondstage_mills = (normpdf(norminv(cdf('logistic',monte_dgp_secondstage_zij,0,1))))./(cdf('logistic',monte_dgp_secondstage_zij,0,1));
zij_real = fixed_effects_selection(monte_dgp_secondstage(:,1)) + fixed_effects_selection(monte_dgp_secondstage(:,2)) + beta_21*monte_dgp_secondstage(:,4) + beta_22*monte_dgp_secondstage(:,5) + beta_23*monte_dgp_secondstage(:,9);
secondstage_mills_real = (normpdf(norminv(cdf('logistic',zij_real,0,1))))./(cdf('logistic',zij_real,0,1));
zij_real_estimated = [monte_dgp_secondstage(:,1),monte_dgp_secondstage(:,2),zij_real,monte_dgp_secondstage_zij];

% run OLS for second stage with FE and including inverse mills ratio
D_i = dummyvar(categorical(monte_dgp_secondstage(:,1)));
D_j = dummyvar(categorical(monte_dgp_secondstage(:,2)));
table_secondstage = array2table([monte_dgp_secondstage(:,4),monte_dgp_secondstage(:,5),monte_dgp_secondstage_mills,D_i,D_j,monte_dgp_secondstage(:,3)]);
second_stage = fitlm(table_secondstage,'Intercept',false);
betas = table2array(second_stage.Coefficients(1:3,1));

table_secondstage_real = array2table([monte_dgp_secondstage(:,4),monte_dgp_secondstage(:,5),secondstage_mills_real,D_i,D_j,monte_dgp_secondstage(:,3)]);
second_stage_real = fitlm(table_secondstage_real,'Intercept',false);
betas_real = table2array(second_stage_real.Coefficients(1:3,1));


function [cmle se m_star combinations ss] = ConditionalLogit(y2,x1,x2,x3)

% construct quadruples
[zz rr1 rr2 rr3 ss combinations] = Rearrangement(y2,x1,x2,x3);  m_star = sum(ss);
% drop non-informative quadruples
zzz = zz(ss==1); zzz = (zzz+1)/2;
rrr1 = rr1(ss==1);
rrr2 = rr2(ss==1);
rrr3 = rr3(ss==1);
rrr = [rrr1 rrr2 rrr3];
% estimate logit on quadruples
[cmle dev stats] = glmfit(rrr,zzz,'binomial','constant','off'); %cmle = coeff(1); se_mle = stats.se(1); t_mle = (mle-theta)./se_mle;
% construct standard errors
[se] = StandardError(cmle,y2,x1,x2,x3,m_star);

function [z r1 r2 r3 s combinations] = Rearrangement(y2,x1,x2,x3)
% auxiliary parameters and functions
F = @(e) 1./(1+exp(-e)); f = @(e) F(e).*(1-F(e)); ff = @(e) f(e).*(1-F(e))-f(e).*F(e); % unit logistic F, f, and f'
[n m] =size(y2); nn = nchoosek(n,2); mm = nchoosek(m-2,2); rho = nn*mm; 

% initialization
z = zeros(rho,1); r1 = z; r2 = z; r3 = z; s = z;

% construction loop (exploiting symmetry in i and j)
c = 1;
for i1 = 1:n,
    for j1=1:m,
        if j1==i1, continue; end
        y2_11 = y2(i1,j1); x1_11 = x1(i1,j1); x2_11 = x2(i1,j1); x3_11 = x3(i1,j1);
        for i2 = i1+1:n,
            if i2==i1 || i2==j1, continue; end
            y2_21 = y2(i2,j1); x1_21 = x1(i2,j1); x2_21 = x2(i2,j1); x3_21 = x3(i2,j1);
            for j2 =  j1+1:m,
                if j2==i1 || j2==j1 || j2==i2, continue; end
                y2_12 = y2(i1,j2); x1_12 = x1(i1,j2); x2_12 = x2(i1,j2); x3_12 = x3(i1,j2);
                y2_22 = y2(i2,j2); x1_22 = x1(i2,j2); x2_22 = x2(i2,j2); x3_22 = x3(i2,j2);

                s(c) = (y2_11>y2_12 & y2_21<y2_22) + (y2_11<y2_12 & y2_21>y2_22);
                z(c) = (y2_11-y2_12)-(y2_21-y2_22)                      ; 
                r1(c) = (x1_11-x1_12)-(x1_21-x1_22)                      ; 
                r2(c) = (x2_11-x2_12)-(x2_21-x2_22)                      ;
                r3(c) = (x3_11-x3_12)-(x3_21-x3_22)                    ;
                
                combinations(c,1) = i1;
                combinations(c,2) = i2;
                combinations(c,3) = j1;
                combinations(c,4) = j2; c = c+1;
            end
        end
    end
end
z = z/2;

function [se] = StandardError(psi,y2,x1,x2,x3,m_star) %changed this function such that it accomodates multivariate
% auxiliary parameters and functions
F = @(e) 1./(1+exp(-e)); f = @(e) F(e).*(1-F(e)); ff = @(e) f(e).*(1-F(e))-f(e).*F(e);  % unit logistic F, f, and f' and f''
[n m] =size(y2); nn = nchoosek(n,2); mm = nchoosek(m-2,2);  rho = n*(n-2)*(n-3)*(n-4); pn = m_star/rho;

dimen = size(psi,1);
S = cell(dimen, 1); %store score for each of the parameters
for d=1:dimen
    S{d} = zeros(n,m);
end
J = zeros(dimen,dimen); %store jacobians
for i1=1:n
    for j1=1:m
        if j1==i1, continue; end
        for i2=1:n
            if i2==i1 || i2==j1, continue; end
            for j2=1:m
                if j2==i2 || j2==i1 || j2==j1, continue; end
                y2_11 = y2(i1,j1);
                x1_11 = x1(i1,j1);
                x2_11 = x2(i1,j1);
                x3_11 = x3(i1,j1);
                y2_21 = y2(i2,j1);
                x1_21 = x1(i2,j1);
                x2_21 = x2(i2,j1);
                x3_21 = x3(i2,j1);
                y2_12 = y2(i1,j2);
                x1_12 = x1(i1,j2);
                x2_12 = x2(i1,j2);                
                x3_12 = x3(i1,j2);
                y2_22 = y2(i2,j2);
                x1_22 = x1(i2,j2);
                x2_22 = x2(i2,j2);                
                x3_22 = x3(i2,j2); 
                
                nondiag = (j1~=i1 & j1~=i2) + (j2~=i1 & j2~=i2);
                c = (y2_11>y2_12 & y2_21<y2_22) + (y2_11<y2_12 & y2_21>y2_22);
                a = (y2_11>y2_12 & y2_21<y2_22)                      ;
                b = (y2_11<y2_12 & y2_21>y2_22)                      ;
                r = zeros(3,1);
                r(1) = (x1_11-x1_12) - (x1_21 - x1_22);
                r(2) = (x2_11-x2_12) - (x2_21 - x2_22);
                r(3) = (x3_11-x3_12) - (x3_21 - x3_22);
                
                %putting items in S
                for d=1:dimen, 
                    S{d}(i1,j1) = S{d}(i1,j1) + 4*r(d) *(a- F(r'*psi))*c*nondiag/((n-2)*(m-3)); 
                end
                
                %putting items in J
                for d=1:dimen,
                    for dd=1:dimen, 
                        J(d,dd) = J(d,dd) + r(d)*r(dd)*(0- f(r'*psi))*c*nondiag/rho; 
                    end
                end
            end
        end
    end
end

xi = cell(dimen,1); V = zeros(dimen,dimen); Q =zeros(dimen,dimen);
xi = S;
for d=1:dimen,
    for dd=1:dimen,
        V(d,dd) = mean(mean(xi{d}.*xi{dd}))              ; % moment asy variance
        Q(d,dd) = J(d,dd); % limit jacobian
    end
end
W = inv(Q'*Q)'*(Q'*V*Q)*inv(Q'*Q); se = sqrt(diag(W)/(n*m)); 