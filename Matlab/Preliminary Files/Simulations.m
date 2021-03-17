%focusing here on design 5

function [OUTPUT]=Simulations
% parameters
R  = 10; % # of Monte Carlo replications
n  =   50; % # of nodes 
beta_11 = 1;  % parameter value
beta_12 = 2.5;
beta_21 = 0.8;
beta_22 = 1;
beta_23 = 2;

% initialization
CMLE = cell(R, 1); SECMLE = CMLE ;%SE_X = CMLE;
% Replications
parfor r=1:R, [CMLE{r} SECMLE{r}] = Simulation(n); end
% prepare screem output
%TCMLE = (CMLE-ones(R,1)*theta)./SECMLE; aCMLE = abs(TCMLE)>=1.96;
%OUTPUT = [CMLE SECMLE TCMLE aCMLE];
OUTPUT = [CMLE SECMLE]

function [cmle se_cmle] = Simulation(n)
warning('off','all'); % prevent overparametrization warning for MLE
% auxiliary parameters and functions
beta_11 = 1;  % parameter value
beta_12 = 2.5;
beta_21 = 0.8;
beta_22 = 1;
beta_23 = 2;
m = nchoosek(n,2)*nchoosek(n-2,2); 
F = @(u) 1./(1+exp(-u)); Q = @(u) log(u./(1-u)); f = @(u) F(u).*(1-F(u)); df = @(u) f(u).*(1-2*F(u));

% data generation % see paper for details on design choices ;

% generating FE
a = repmat(1:(n*n),1)';
monte_fe_i_index = sort(a-n*floor(a/n) + 1); 
monte_fe_j_index = (a-n*floor(a/n) + n + 1); 
monte_fe_indices = [monte_fe_i_index,monte_fe_j_index];
monte_fe_indices = sortrows(monte_fe_indices, [1 2]);
reference = monte_fe_indices(:, 2)-n;
monte_fe_indices = [monte_fe_indices,reference];
condition = monte_fe_indices(:,1)==monte_fe_indices(:,3);
monte_fe_indices(condition,:) = [];
monte_fe_indices(:,3)=[];
fixed_effects_vartheta_i = zeros(size(monte_fe_indices,1),1);
fixed_effects_chi_j = zeros(size(monte_fe_indices,1),1);
fixed_effects_xi_i = zeros(size(monte_fe_indices,1),1);
fixed_effects_zeta_j = zeros(size(monte_fe_indices,1),1);
monte_fe_indices = [monte_fe_indices,fixed_effects_vartheta_i, fixed_effects_chi_j, fixed_effects_xi_i,fixed_effects_zeta_j];
fixed_effects_observation = normrnd(0, 1, [(n+n), 1]);
fixed_effects_selection = normrnd(0, 1, [(n+n), 1]);
for i=1:size(monte_fe_indices,1)
    monte_fe_indices(i,3)= fixed_effects_observation(monte_fe_indices(i,1));
    monte_fe_indices(i,4)= fixed_effects_observation(monte_fe_indices(i,2));
    monte_fe_indices(i,5)= fixed_effects_selection(monte_fe_indices(i,1));
    monte_fe_indices(i,6)= fixed_effects_selection(monte_fe_indices(i,2));
end
monte_fe_indices(:,3) = monte_fe_indices(:,3) - mean(monte_fe_indices(:,3));
monte_fe_indices(:,4) = monte_fe_indices(:,4) - mean(monte_fe_indices(:,4));
monte_fe_indices(:,5) = monte_fe_indices(:,5) - mean(monte_fe_indices(:,5));
monte_fe_indices(:,6) = monte_fe_indices(:,6) - mean(monte_fe_indices(:,6));

% generating explanatory variables
number_observations = n*(n-1);
monte_x1 = normrnd(0, 1, [number_observations, 1]) + monte_fe_indices(:,5) + monte_fe_indices(:,6);
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
monte_y2 = (beta_21*monte_x1 + beta_22*monte_x2  + beta_23*monte_x3 + monte_fe_indices(:,5) + monte_fe_indices(:,6) + errors_selection)>=0;
monte_y1_x = beta_11*monte_x1 + beta_12*monte_x2 + monte_fe_indices(:,3) + monte_fe_indices(:,4) +  errors_observation;
monte_y1 = monte_y1_x.*(monte_y2>0);
monte_y1(monte_y1==0) = NaN;

%put all in a matrix
monte_dgp = [monte_fe_indices(:,1),monte_fe_indices(:,2),monte_y1,monte_fe_indices(:,3),monte_fe_indices(:,4),monte_x1, monte_x2,monte_y2,monte_fe_indices(:,5),monte_fe_indices(:,6),monte_x3];

%create matrices for conditional logit estimation
%Put back diagonal elements
add =  NaN(n,size(monte_dgp,2));
add1 = repmat(1:n,1)';
add2 = add1+n;
add(:,1)= add1;
add(:,2)= add2;
monte_dgp = vertcat(monte_dgp,add);
%create ij variable
monte_dgp = sortrows(monte_dgp, [1 2]);
ij_monte_dgp = strcat(num2str(monte_dgp(:,1)),num2str(monte_dgp(:,2)));
ij_monte_dgp = mat2cell(ij_monte_dgp, ones(size(ij_monte_dgp,1),1), size(ij_monte_dgp,2));
ij_monte_dgp = cellfun(@str2num, ij_monte_dgp,'UniformOutput',false);
%create matrices
y_matrix_check_indices = reshape(ij_monte_dgp,n,n)';
monte_y2_matrix = reshape(monte_dgp(:,8),n,n)';
monte_x1_matrix = reshape(monte_dgp(:,6),n,n)';
monte_x2_matrix = reshape(monte_dgp(:,7),n,n)';
monte_x3_matrix = reshape(monte_dgp(:,11),n,n)';

% conditional logit estimation
[cmle se_cmle m_star] = ConditionalLogit(monte_y2_matrix,monte_x1_matrix,monte_x2_matrix,monte_x3_matrix); 
%check if in original file if there is something after this


function [cmle se m_star] = ConditionalLogit(y2,x1,x2,x3)

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
[se] = StandardError(cmle,y2,x1,x2,x3,m_star);

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
                r3(c) = (x3_11-x3_12)-(x3_21-x3_22)                      ; c = c+1;
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




