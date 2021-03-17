%%%% Author: Gabriela Szini
%%%% Date: 28/07/2020
%%%% File based on code in
%%%% Prelimiary/Simulations_Jochmans_Hybrid_correctFE_optimize

%%%%%%%%%%%%%%%%%%%%%%%%% DESIGN 9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [OUTPUT]=Simulations
% parameters
R  = 500; % # of Monte Carlo replications
n  =  25; % # of nodes 
beta_11 = 1;  % parameter value
beta_12 = 2.5;
beta_21 = 0.8;
beta_22 = 1;
beta_23 = 2;

% initialization
CMLE = cell(R, 1); SECMLE = CMLE ; BETAS = CMLE;  FE_SELECTION = CMLE; ZIJ = CMLE;
% Replications
parfor r=1:R, [CMLE{r} SECMLE{r} BETAS{r} FE_SELECTION{r} ZIJ{r}] = Simulation(n, beta_11, beta_12, beta_21, beta_22, beta_23); end
% prepare screem output
%TCMLE = (CMLE-ones(R,1)*theta)./SECMLE; aCMLE = abs(TCMLE)>=1.96;
%OUTPUT = [CMLE SECMLE TCMLE aCMLE];
OUTPUT = [CMLE SECMLE BETAS FE_SELECTION ZIJ]
save('Simulations_Jochmans_Hybrid_Design9_n25.mat','OUTPUT')

function [cmle se_cmle betas FE_selection_real_estimated zij_real_estimated] = Simulation(n, beta_11, beta_12, beta_21, beta_22, beta_23)
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
monte_fe_indices = [monte_fe_indices,zeros(size(monte_fe_indices,1),4)]; %here as well
fixed_effects_observation = normrnd(0, 1, [(n+n), 1]);
fixed_effects_selection = normrnd(0, 1, [(n+n), 1]);
monte_fe_indices(:,3)= fixed_effects_observation(monte_fe_indices(:,1)) - mean(fixed_effects_observation(monte_fe_indices(:,1)));
monte_fe_indices(:,4)= fixed_effects_observation(monte_fe_indices(:,2)) - mean(fixed_effects_observation(monte_fe_indices(:,2)));
monte_fe_indices(:,5)= fixed_effects_selection(monte_fe_indices(:,1)) - mean(fixed_effects_selection(monte_fe_indices(:,1)));
monte_fe_indices(:,6)= fixed_effects_selection(monte_fe_indices(:,2)) - mean(fixed_effects_selection(monte_fe_indices(:,2)));

% generating explanatory variables
number_observations = n*(n-1);
monte_x1 = normrnd(0, 1, [number_observations, 1]) + monte_fe_indices(:,5) + monte_fe_indices(:,6) + monte_fe_indices(:,3) + monte_fe_indices(:,4);
monte_x2 = (rand(number_observations,1)) <=0.5;
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
monte_y2 = (beta_21*monte_x1 + beta_22*monte_x2  + beta_23*monte_x3 + monte_fe_indices(:,5) + monte_fe_indices(:,6) + errors_selection)>=0;
monte_y1 = (beta_11*monte_x1 + beta_12*monte_x2 + monte_fe_indices(:,3) + monte_fe_indices(:,4) +  errors_observation).*(monte_y2>0);
monte_y1(monte_y1==0) = NaN;

%put all in a matrix
monte_dgp = [monte_fe_indices(:,1),monte_fe_indices(:,2),monte_y1,monte_fe_indices(:,3),monte_fe_indices(:,4),monte_x1, monte_x2,monte_y2,monte_fe_indices(:,5),monte_fe_indices(:,6),monte_x3];

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
monte_y2_matrix = reshape(monte_dgp(:,8),n,n)';
monte_x1_matrix = reshape(monte_dgp(:,6),n,n)';
monte_x2_matrix = reshape(monte_dgp(:,7),n,n)';
monte_x3_matrix = reshape(monte_dgp(:,11),n,n)';
monte_y1_matrix = reshape(monte_dgp(:,3),n,n)';

% conditional logit estimation - changed here as we need more outputs
% coming from se
[cmle se_cmle m_star x1_k x2_k x3_k y1_k y2_k combinations] = ConditionalLogit(monte_y2_matrix,monte_x1_matrix,monte_x2_matrix,monte_x3_matrix,monte_y1_matrix); 

%%%% Kyriazidou estimation %%%%

% turning transformed data into column vectors
x1_k = x1_k'; x2_k = x2_k'; x3_k = x3_k'; y1_k = y1_k'; y2_k = y2_k';

% defining parameters first set of parameters, h=1
r=1; h=1; delta=0.1;
% h = 0.5, 2, 3

% for true second stage we set
cmle_1 = beta_21; cmle_2 = beta_22; cmle_3 = beta_23;

[h_nc_true_1 h_nc_delta_true_1 betas_true_1 betas_corrected_true_1] = Kyriazidou(cmle_1, cmle_2, cmle_3, x1_k, x2_k, x3_k, y1_k, y2_k, combinations, r, h, delta);
%then do the same function for estimated parameters and other h's
%check different nc : only non-missing


function [h_nc h_nc_delta betas betas_corrected] = Kyriazidou(cmle_1, cmle_2, cmle_3, x1_k, x2_k, x3_k, y1_k, y2_k, combinations, r, h, delta)
nc = size(combinations,1);
h_nc = h*(nc^(-1/(2*(r+1)+1)));
h_nc_delta = h*(nc^(-delta/(2*(r+1)+1)));
% construct Y1_ij, X1_ij and X2_ij
weights = sqrt(normpdf((cmle_1*x1_k + cmle_2*x2_k + cmle_3*x3_k)/h_nc));
Y1_ij = weights.*y1_k;
X1_ij = weights.*x1_k;
X2_ij = weights.*x2_k;
table_secondstage = array2table([X1_ij,X2_ij,Y1_ij]);
second_stage = fitlm(table_secondstage,'Intercept',false);
betas = table2array(second_stage.Coefficients(1:2,1));
% construct Y1_ij_delta, X1_ij_delta and X2_ij_delta
weights_delta = sqrt(normpdf((cmle_1*x1_k + cmle_2*x2_k + cmle_3*x3_k)/h_nc_delta));
Y1_ij_delta = weights_delta.*y1_k;
X1_ij_delta = weights_delta.*x1_k;
X2_ij_delta = weights_delta.*x2_k;
table_secondstage_delta = array2table([X1_ij_delta,X2_ij_delta,Y1_ij_delta]);
second_stage_delta = fitlm(table_secondstage_delta,'Intercept',false);
betas_delta = table2array(second_stage_delta.Coefficients(1:2,1));
betas_corrected = (betas - nc^(-(1-delta)*(r+1)/(2*(r+1)+1))*betas_delta)/(1 - nc^(-(1-delta)*(r+1)/(2*(r+1)+1)));

function [cmle se m_star x1_k x2_k x3_k y1_k y2_k combinations] = ConditionalLogit(y2,x1,x2,x3,y1)

% construct quadruples
[zz rr1 rr2 rr3 ss] = Rearrangement(y2,x1,x2,x3);  m_star = sum(ss);
% drop non-informative quadruples
zzz = zz(ss==1); zzz = (zzz+1)/2;
rrr1 = rr1(ss==1);
rrr2 = rr2(ss==1);
rrr3 = rr3(ss==1);
rrr = [rrr1 rrr2 rrr3];
% estimate logit on quadruples
[cmle dev stats] = glmfit(rrr,zzz,'binomial','constant','off'); %cmle = coeff(1); se_mle = stats.se(1); t_mle = (mle-theta)./se_mle;
% construct standard errors
% we take all possible combinations from the se function as in this way do
% not take into account symmetry in i and j
[se x1_k x2_k x3_k y1_k y2_k combinations] = StandardError(cmle,y2,x1,x2,x3,m_star,y1);

function [z r1 r2 r3 s] = Rearrangement(y2,x1,x2,x3)
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
                
                combinations_jochmans(c,1) = i1;
                combinations_jochmans(c,2) = i2;
                combinations_jochmans(c,3) = j1;
                combinations_jochmans(c,4) = j2; c = c+1;
            end
        end
    end
end
z = z/2;

function [se x1_k x2_k x3_k y1_k y2_k combinations] = StandardError(psi,y2,x1,x2,x3,m_star,y1) %changed this function such that it accomodates multivariate
% auxiliary parameters and functions
F = @(e) 1./(1+exp(-e)); f = @(e) F(e).*(1-F(e)); ff = @(e) f(e).*(1-F(e))-f(e).*F(e);  % unit logistic F, f, and f' and f''
[n m] =size(y2); nn = nchoosek(n,2); mm = nchoosek(m-2,2);  rho = n*(n-2)*(n-3)*(n-4); pn = m_star/rho;

dimen = size(psi,1);
S = cell(dimen, 1); %store score for each of the parameters
for d=1:dimen
    S{d} = zeros(n,m);
end
J = zeros(dimen,dimen); %store jacobians
count = 1;
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
                y1_11 = y1(i1,j1);
                y2_21 = y2(i2,j1);
                x1_21 = x1(i2,j1);
                x2_21 = x2(i2,j1);
                x3_21 = x3(i2,j1);
                y1_21 = y1(i2,j1);
                y2_12 = y2(i1,j2);
                x1_12 = x1(i1,j2);
                x2_12 = x2(i1,j2);                
                x3_12 = x3(i1,j2);
                y1_12 = y1(i1,j2);
                y2_22 = y2(i2,j2);
                x1_22 = x1(i2,j2);
                x2_22 = x2(i2,j2);                
                x3_22 = x3(i2,j2); 
                y1_22 = y1(i2,j2); 
                
                nondiag = (j1~=i1 & j1~=i2) + (j2~=i1 & j2~=i2);
                c = (y2_11>y2_12 & y2_21<y2_22) + (y2_11<y2_12 & y2_21>y2_22);
                a = (y2_11>y2_12 & y2_21<y2_22)                      ;
                b = (y2_11<y2_12 & y2_21>y2_22)                      ;
                r = zeros(3,1);
                r(1) = (x1_11-x1_12) - (x1_21 - x1_22);
                r(2) = (x2_11-x2_12) - (x2_21 - x2_22);
                r(3) = (x3_11-x3_12) - (x3_21 - x3_22);
                
                %storing for kyriazidou approach:
                x1_k(count) = r(1);
                x2_k(count) = r(2);
                x3_k(count) = r(3);
                y1_k(count) = (y1_11-y1_12) - (y1_21 - y1_22);
                y2_k(count) = (y2_11-y2_12) - (y2_21 - y2_22);
                combinations(count,1) = i1;
                combinations(count,2) = i2;
                combinations(count,3) = j1;
                combinations(count,4) = j2; count = count+1;
                
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
