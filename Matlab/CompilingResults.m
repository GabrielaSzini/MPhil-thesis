%%%% Author: Gabriela Szini
%%%% Date: 01/08/2020

%%%%%%%%%%%%%%%%%%%%%% COMPILING RESULTS FOR TABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% Hybrid approach %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Design 1 / 3 / 6
clearvars 
load('Simulations_Jochmans_Hybrid_Design6_n25.mat')

betas_first_stage_hybrid = cell2mat(OUTPUT(:,1)');
save('/Users/gabrielaszini/Documents/Master Thesis/Simulations/Figures/betas_first_stage_charbonneau_Design6.mat','betas_first_stage_hybrid')
% bias
beta1_bias_first_stage_hybrid = mean(betas_first_stage_hybrid(1,:))-0.8
beta2_bias_first_stage_hybrid = mean(betas_first_stage_hybrid(2,:))-1
beta3_bias_first_stage_hybrid = mean(betas_first_stage_hybrid(3,:))-2
% standard deviation
beta1_std_first_stage_hybrid = std(betas_first_stage_hybrid(1,:))
beta2_std_first_stage_hybrid = std(betas_first_stage_hybrid(2,:))
beta3_std_first_stage_hybrid = std(betas_first_stage_hybrid(3,:))
% mean standard error
se_first_stage_hybrid = cell2mat(OUTPUT(:,2)');
beta1_se_first_stage_hybrid = mean(se_first_stage_hybrid(1,:))
beta2_se_first_stage_hybrid = mean(se_first_stage_hybrid(2,:))
beta3_se_first_stage_hybrid = mean(se_first_stage_hybrid(3,:))
save('/Users/gabrielaszini/Documents/Master Thesis/Simulations/Figures/se_first_stage_charbonneau_Design6.mat','se_first_stage_hybrid')
% test size
ttest_beta1_first_stage_hybrid = (betas_first_stage_hybrid(1,:)-0.8)./se_first_stage_hybrid(1,:);
size_ttest_beta1_first_stage_hybrid = sum(abs(ttest_beta1_first_stage_hybrid)>1.96)/size(betas_first_stage_hybrid,2)
ttest_beta2_first_stage_hybrid = (betas_first_stage_hybrid(2,:)-1)./se_first_stage_hybrid(2,:);
size_ttest_beta2_first_stage_hybrid = sum(abs(ttest_beta2_first_stage_hybrid)>1.96)/size(betas_first_stage_hybrid,2)
ttest_beta3_first_stage_hybrid = (betas_first_stage_hybrid(3,:)-2)./se_first_stage_hybrid(3,:);
size_ttest_beta3_first_stage_hybrid = sum(abs(ttest_beta3_first_stage_hybrid)>1.96)/size(betas_first_stage_hybrid,2)

betas_second_stage_hybrid = cell2mat(OUTPUT(:,3)');
% bias 
beta1_bias_second_stage_hybrid = mean(betas_second_stage_hybrid(1,:))-1
beta2_bias_second_stage_hybrid = mean(betas_second_stage_hybrid(2,:))-2.5
beta3_bias_second_stage_hybrid = mean(betas_second_stage_hybrid(3,:))-(-0.7)
save('/Users/gabrielaszini/Documents/Master Thesis/Simulations/Figures/betas_second_stage_hybrid_Design6.mat','betas_second_stage_hybrid')

% standard deviation
beta1_std_second_stage_hybrid = std(betas_second_stage_hybrid(1,:))
beta2_std_second_stage_hybrid = std(betas_second_stage_hybrid(2,:))
beta3_std_second_stage_hybrid = std(betas_second_stage_hybrid(3,:))


% zij : first real then estimated
zij_real_estimated_hybrid = cell2mat(OUTPUT(:,4));
z12_real_estimated_hybrid = zij_real_estimated_hybrid;
condition1 = z12_real_estimated_hybrid(:,1) ~= 1;
z12_real_estimated_hybrid(condition1,:)=[];  
condition2 = z12_real_estimated_hybrid(:,2) ~= 27;
z12_real_estimated_hybrid(condition2,:)=[];  
z12_bias_hybrid = mean(z12_real_estimated_hybrid(:,4)-z12_real_estimated_hybrid(:,3))
z12_std_hybrid = std(z12_real_estimated_hybrid(:,4)-z12_real_estimated_hybrid(:,3)) %%%UPDATED HERE
save('/Users/gabrielaszini/Documents/Master Thesis/Simulations/Figures/zij_hybrid_Design6.mat','z12_real_estimated_hybrid')

betas_real_second_stage_hybrid = cell2mat(OUTPUT(:,5)');
% bias 
beta1_real_bias_second_stage_hybrid = mean(betas_real_second_stage_hybrid(1,:))-1
beta2_real_bias_second_stage_hybrid = mean(betas_real_second_stage_hybrid(2,:))-2.5
beta3_real_bias_second_stage_hybrid = mean(betas_real_second_stage_hybrid(3,:))-(-0.7)

% standard deviation
beta1_real_std_second_stage_hybrid = std(betas_real_second_stage_hybrid(1,:))
beta2_real_std_second_stage_hybrid = std(betas_real_second_stage_hybrid(2,:))
beta3_real_std_second_stage_hybrid = std(betas_real_second_stage_hybrid(3,:))


clearvars
%%% Design 2 /4 /5 /9
load('Simulations_Jochmans_Hybrid_Design9_n25.mat')

betas_first_stage_hybrid = cell2mat(OUTPUT(:,1)');
% bias
beta1_bias_first_stage_hybrid = mean(betas_first_stage_hybrid(1,:))-0.8
beta2_bias_first_stage_hybrid = mean(betas_first_stage_hybrid(2,:))-1
beta3_bias_first_stage_hybrid = mean(betas_first_stage_hybrid(3,:))-2
% standard deviation
beta1_std_first_stage_hybrid = std(betas_first_stage_hybrid(1,:))
beta2_std_first_stage_hybrid = std(betas_first_stage_hybrid(2,:))
beta3_std_first_stage_hybrid = std(betas_first_stage_hybrid(3,:))
save('/Users/gabrielaszini/Documents/Master Thesis/Simulations/Figures/betas_first_stage_charbonneau_Design7.mat','betas_first_stage_hybrid')

% mean standard error
se_first_stage_hybrid = cell2mat(OUTPUT(:,2)');
beta1_se_first_stage_hybrid = mean(se_first_stage_hybrid(1,:))
beta2_se_first_stage_hybrid = mean(se_first_stage_hybrid(2,:))
beta3_se_first_stage_hybrid = mean(se_first_stage_hybrid(3,:))
save('/Users/gabrielaszini/Documents/Master Thesis/Simulations/Figures/se_first_stage_charbonneau_Design7.mat','se_first_stage_hybrid')

% test size
ttest_beta1_first_stage_hybrid = (betas_first_stage_hybrid(1,:)-0.8)./se_first_stage_hybrid(1,:);
size_ttest_beta1_first_stage_hybrid = sum(abs(ttest_beta1_first_stage_hybrid)>1.96)/size(betas_first_stage_hybrid,2)
ttest_beta2_first_stage_hybrid = (betas_first_stage_hybrid(2,:)-1)./se_first_stage_hybrid(2,:);
size_ttest_beta2_first_stage_hybrid = sum(abs(ttest_beta2_first_stage_hybrid)>1.96)/size(betas_first_stage_hybrid,2)
ttest_beta3_first_stage_hybrid = (betas_first_stage_hybrid(3,:)-2)./se_first_stage_hybrid(3,:);
size_ttest_beta3_first_stage_hybrid = sum(abs(ttest_beta3_first_stage_hybrid)>1.96)/size(betas_first_stage_hybrid,2)

betas_second_stage_hybrid = cell2mat(OUTPUT(:,3)');
% bias 
beta1_bias_second_stage_hybrid = mean(betas_second_stage_hybrid(1,:))-1
beta2_bias_second_stage_hybrid = mean(betas_second_stage_hybrid(2,:))-2.5
beta3_bias_second_stage_hybrid = mean(betas_second_stage_hybrid(3,:))-(-0.7)
save('/Users/gabrielaszini/Documents/Master Thesis/Simulations/Figures/betas_second_stage_hybrid_Design7.mat','betas_second_stage_hybrid')

% standard deviation
beta1_std_second_stage_hybrid = std(betas_second_stage_hybrid(1,:))
beta2_std_second_stage_hybrid = std(betas_second_stage_hybrid(2,:))
beta3_std_second_stage_hybrid = std(betas_second_stage_hybrid(3,:))


% fixed effects
fe_real_estimated_hybrid = cell2mat(OUTPUT(:,4));
fe_1_real_estimated_hybrid = fe_real_estimated_hybrid;
fe_1_real_estimated_hybrid(fe_1_real_estimated_hybrid(:,1) ~= 1,:)=[];  
fe_1_bias = nanmean(fe_1_real_estimated_hybrid(:,3)-fe_1_real_estimated_hybrid(:,2))
fe_1_std = nanstd(fe_1_real_estimated_hybrid(:,3)-fe_1_real_estimated_hybrid(:,2))
fe_2_real_estimated_hybrid = fe_real_estimated_hybrid;
fe_2_real_estimated_hybrid(fe_2_real_estimated_hybrid(:,1) ~= 27,:)=[];  
fe_2_bias = nanmean(fe_2_real_estimated_hybrid(:,3)-fe_2_real_estimated_hybrid(:,2))
fe_2_std = nanstd(fe_2_real_estimated_hybrid(:,3)-fe_2_real_estimated_hybrid(:,2))
save('/Users/gabrielaszini/Documents/Master Thesis/Simulations/Figures/fe1_hybrid_Design7.mat','fe_1_real_estimated_hybrid')
save('/Users/gabrielaszini/Documents/Master Thesis/Simulations/Figures/fe2_hybrid_Design7.mat','fe_2_real_estimated_hybrid')


% zij : first real then estimated
zij_real_estimated_hybrid = cell2mat(OUTPUT(:,5));
z12_real_estimated_hybrid = zij_real_estimated_hybrid;
condition1 = z12_real_estimated_hybrid(:,1) ~= 1;
z12_real_estimated_hybrid(condition1,:)=[];  
condition2 = z12_real_estimated_hybrid(:,2) ~= 27;
z12_real_estimated_hybrid(condition2,:)=[];  
z12_bias_hybrid = mean(z12_real_estimated_hybrid(:,4)-z12_real_estimated_hybrid(:,3))
z12_std = std(z12_real_estimated_hybrid(:,4)-z12_real_estimated_hybrid(:,3))
save('/Users/gabrielaszini/Documents/Master Thesis/Simulations/Figures/zij_hybrid_Design7.mat','z12_real_estimated_hybrid')

betas_real_second_stage_hybrid = cell2mat(OUTPUT(:,6)');
% bias 
beta1_real_bias_second_stage_hybrid = mean(betas_real_second_stage_hybrid(1,:))-1
beta2_real_bias_second_stage_hybrid = mean(betas_real_second_stage_hybrid(2,:))-2.5
beta3_real_bias_second_stage_hybrid = mean(betas_real_second_stage_hybrid(3,:))-(-0.7)

% standard deviation
beta1_real_std_second_stage_hybrid = std(betas_real_second_stage_hybrid(1,:))
beta2_real_std_second_stage_hybrid = std(betas_real_second_stage_hybrid(2,:))
beta3_real_std_second_stage_hybrid = std(betas_real_second_stage_hybrid(3,:))


%%%%%%%%%%%%%%%%%%%%%%%% Kyriazidou approach %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars 
load('Simulations_Jochmans_Kyriazidou_Design5_n25.mat')

betas_first_stage_kyri = cell2mat(OUTPUT(:,1)');
% bias
beta1_bias_first_stage_kyri = mean(betas_first_stage_kyri(1,:))-0.8
beta2_bias_first_stage_kyri = mean(betas_first_stage_kyri(2,:))-1
beta3_bias_first_stage_kyri = mean(betas_first_stage_kyri(3,:))-2
% standard deviation
beta1_std_first_stage_kyri = std(betas_first_stage_kyri(1,:))
beta2_std_first_stage_kyri = std(betas_first_stage_kyri(2,:))
beta3_std_first_stage_kyri = std(betas_first_stage_kyri(3,:))
% mean standard error
se_first_stage_kyri = cell2mat(OUTPUT(:,2)');
beta1_se_first_stage_kyri = mean(se_first_stage_kyri(1,:))
beta2_se_first_stage_kyri = mean(se_first_stage_kyri(2,:))
beta3_se_first_stage_kyri = mean(se_first_stage_kyri(3,:))
% test size
ttest_beta1_first_stage_kyri = (betas_first_stage_kyri(1,:)-0.8)./se_first_stage_kyri(1,:);
size_ttest_beta1_first_stage_kyri = sum(abs(ttest_beta1_first_stage_kyri)>1.96)/size(betas_first_stage_kyri,2)
ttest_beta2_first_stage_kyri = (betas_first_stage_kyri(2,:)-1)./se_first_stage_kyri(2,:);
size_ttest_beta2_first_stage_kyri = sum(abs(ttest_beta2_first_stage_kyri)>1.96)/size(betas_first_stage_kyri,2)
ttest_beta3_first_stage_kyri = (betas_first_stage_kyri(3,:)-2)./se_first_stage_kyri(3,:);
size_ttest_beta3_first_stage_kyri = sum(abs(ttest_beta3_first_stage_kyri)>1.96)/size(betas_first_stage_kyri,2)

% with true first stage
betas_second_stage_true_kyri = cell2mat(OUTPUT(:,5));
betas_second_stage_true_kyri_05 = betas_second_stage_true_kyri(:,1);
betas_second_stage_true_kyri_05 = reshape(betas_second_stage_true_kyri_05, 2,500);
beta1_second_stage_true_kyri_05 = mean(betas_second_stage_true_kyri_05(1,:))-1
beta1_std_second_stage_true_kyri_05 = std(betas_second_stage_true_kyri_05(1,:))
beta2_second_stage_true_kyri_05 = mean(betas_second_stage_true_kyri_05(2,:))-2.5
beta2_std_second_stage_true_kyri_05 = std(betas_second_stage_true_kyri_05(2,:))

betas_second_stage_true_kyri_1 = betas_second_stage_true_kyri(:,2);
betas_second_stage_true_kyri_1 = reshape(betas_second_stage_true_kyri_1, 2,500);
beta1_second_stage_true_kyri_1 = mean(betas_second_stage_true_kyri_1(1,:))-1
beta1_std_second_stage_true_kyri_1 = std(betas_second_stage_true_kyri_1(1,:))
beta2_second_stage_true_kyri_1 = mean(betas_second_stage_true_kyri_1(2,:))-2.5
beta2_std_second_stage_true_kyri_1 = std(betas_second_stage_true_kyri_1(2,:))

betas_second_stage_true_kyri_2 = betas_second_stage_true_kyri(:,3);
betas_second_stage_true_kyri_2 = reshape(betas_second_stage_true_kyri_2, 2,500);
beta1_second_stage_true_kyri_2 = mean(betas_second_stage_true_kyri_2(1,:))-1
beta1_std_second_stage_true_kyri_2 = std(betas_second_stage_true_kyri_2(1,:))
beta2_second_stage_true_kyri_2 = mean(betas_second_stage_true_kyri_2(2,:))-2.5
beta2_std_second_stage_true_kyri_2 = std(betas_second_stage_true_kyri_2(2,:))

betas_second_stage_true_kyri_3 = betas_second_stage_true_kyri(:,4);
betas_second_stage_true_kyri_3 = reshape(betas_second_stage_true_kyri_3, 2,500);
beta1_second_stage_true_kyri_3 = mean(betas_second_stage_true_kyri_3(1,:))-1
beta1_std_second_stage_true_kyri_3 = std(betas_second_stage_true_kyri_3(1,:))
beta2_second_stage_true_kyri_3 = mean(betas_second_stage_true_kyri_3(2,:))-2.5
beta2_std_second_stage_true_kyri_3 = std(betas_second_stage_true_kyri_3(2,:))

betas_second_stage_true_kyri_5 = betas_second_stage_true_kyri(:,5);
betas_second_stage_true_kyri_5 = reshape(betas_second_stage_true_kyri_5, 2,500);
beta1_second_stage_true_kyri_5 = mean(betas_second_stage_true_kyri_5(1,:))-1
beta1_std_second_stage_true_kyri_5 = std(betas_second_stage_true_kyri_5(1,:))
beta2_second_stage_true_kyri_5 = mean(betas_second_stage_true_kyri_5(2,:))-2.5
beta2_std_second_stage_true_kyri_5 = std(betas_second_stage_true_kyri_5(2,:))

betas_second_stage_true_kyri_10 = betas_second_stage_true_kyri(:,6);
betas_second_stage_true_kyri_10 = reshape(betas_second_stage_true_kyri_10, 2,500);
beta1_second_stage_true_kyri_10 = mean(betas_second_stage_true_kyri_10(1,:))-1
beta1_std_second_stage_true_kyri_10 = std(betas_second_stage_true_kyri_10(1,:))
beta2_second_stage_true_kyri_10 = mean(betas_second_stage_true_kyri_10(2,:))-2.5
beta2_std_second_stage_true_kyri_10 = std(betas_second_stage_true_kyri_10(2,:))

% with true first stage and corrected
betas_second_stage_true_corrected_kyri = cell2mat(OUTPUT(:,6));
betas_second_stage_true_corrected_kyri_05 = betas_second_stage_true_corrected_kyri(:,1);
betas_second_stage_true_corrected_kyri_05 = reshape(betas_second_stage_true_corrected_kyri_05, 2,500);
beta1_second_stage_true_corrected_kyri_05 = mean(betas_second_stage_true_corrected_kyri_05(1,:))-1
beta1_std_second_stage_true_corrected_kyri_05 = std(betas_second_stage_true_corrected_kyri_05(1,:))
beta2_second_stage_true_corrected_kyri_05 = mean(betas_second_stage_true_corrected_kyri_05(2,:))-2.5
beta2_std_second_stage_true_corrected_kyri_05 = std(betas_second_stage_true_corrected_kyri_05(2,:))

betas_second_stage_true_corrected_kyri_1 = betas_second_stage_true_corrected_kyri(:,2);
betas_second_stage_true_corrected_kyri_1 = reshape(betas_second_stage_true_corrected_kyri_1, 2,500);
beta1_second_stage_true_corrected_kyri_1 = mean(betas_second_stage_true_corrected_kyri_1(1,:))-1
beta1_std_second_stage_true_corrected_kyri_1 = std(betas_second_stage_true_corrected_kyri_1(1,:))
beta2_second_stage_true_corrected_kyri_1 = mean(betas_second_stage_true_corrected_kyri_1(2,:))-2.5
beta2_std_second_stage_true_corrected_kyri_1 = std(betas_second_stage_true_corrected_kyri_1(2,:))

betas_second_stage_true_corrected_kyri_2 = betas_second_stage_true_corrected_kyri(:,3);
betas_second_stage_true_corrected_kyri_2 = reshape(betas_second_stage_true_corrected_kyri_2, 2,500);
beta1_second_stage_true_corrected_kyri_2 = mean(betas_second_stage_true_corrected_kyri_2(1,:))-1
beta1_std_second_stage_true_corrected_kyri_2 = std(betas_second_stage_true_corrected_kyri_2(1,:))
beta2_second_stage_true_corrected_kyri_2 = mean(betas_second_stage_true_corrected_kyri_2(2,:))-2.5
beta2_std_second_stage_true_corrected_kyri_2 = std(betas_second_stage_true_corrected_kyri_2(2,:))

betas_second_stage_true_corrected_kyri_3 = betas_second_stage_true_corrected_kyri(:,4);
betas_second_stage_true_corrected_kyri_3 = reshape(betas_second_stage_true_corrected_kyri_3, 2,500);
beta1_second_stage_true_corrected_kyri_3 = mean(betas_second_stage_true_corrected_kyri_3(1,:))-1
beta1_std_second_stage_true_corrected_kyri_3 = std(betas_second_stage_true_corrected_kyri_3(1,:))
beta2_second_stage_true_corrected_kyri_3 = mean(betas_second_stage_true_corrected_kyri_3(2,:))-2.5
beta2_std_second_stage_true_corrected_kyri_3 = std(betas_second_stage_true_corrected_kyri_3(2,:))

betas_second_stage_true_corrected_kyri_5 = betas_second_stage_true_corrected_kyri(:,5);
betas_second_stage_true_corrected_kyri_5 = reshape(betas_second_stage_true_corrected_kyri_5, 2,500);
beta1_second_stage_true_corrected_kyri_5 = mean(betas_second_stage_true_corrected_kyri_5(1,:))-1
beta1_std_second_stage_true_corrected_kyri_5 = std(betas_second_stage_true_corrected_kyri_5(1,:))
beta2_second_stage_true_corrected_kyri_5 = mean(betas_second_stage_true_corrected_kyri_5(2,:))-2.5
beta2_std_second_stage_true_corrected_kyri_5 = std(betas_second_stage_true_corrected_kyri_5(2,:))

betas_second_stage_true_corrected_kyri_10 = betas_second_stage_true_corrected_kyri(:,6);
betas_second_stage_true_corrected_kyri_10 = reshape(betas_second_stage_true_corrected_kyri_10, 2,500);
beta1_second_stage_true_corrected_kyri_10 = mean(betas_second_stage_true_corrected_kyri_10(1,:))-1
beta1_std_second_stage_true_corrected_kyri_10 = std(betas_second_stage_true_corrected_kyri_10(1,:))
beta2_second_stage_true_corrected_kyri_10 = mean(betas_second_stage_true_corrected_kyri_10(2,:))-2.5
beta2_std_second_stage_true_corrected_kyri_10 = std(betas_second_stage_true_corrected_kyri_10(2,:))

% with estimated first stage
betas_second_stage_estimated_kyri = cell2mat(OUTPUT(:,9));
betas_second_stage_estimated_kyri_05 = betas_second_stage_estimated_kyri(:,1);
betas_second_stage_estimated_kyri_05 = reshape(betas_second_stage_estimated_kyri_05, 2,500);
save('/Users/gabrielaszini/Documents/Master Thesis/Simulations/Figures/beta_05_kyriazidou_Design5.mat','betas_second_stage_estimated_kyri_05')
beta1_second_stage_estimated_kyri_05 = mean(betas_second_stage_estimated_kyri_05(1,:))-1
beta1_std_second_stage_estimated_kyri_05 = std(betas_second_stage_estimated_kyri_05(1,:))
beta2_second_stage_estimated_kyri_05 = mean(betas_second_stage_estimated_kyri_05(2,:))-2.5
beta2_std_second_stage_estimated_kyri_05 = std(betas_second_stage_estimated_kyri_05(2,:))

betas_second_stage_estimated_kyri_1 = betas_second_stage_estimated_kyri(:,2);
betas_second_stage_estimated_kyri_1 = reshape(betas_second_stage_estimated_kyri_1, 2,500);
beta1_second_stage_estimated_kyri_1 = mean(betas_second_stage_estimated_kyri_1(1,:))-1
beta1_std_second_stage_estimated_kyri_1 = std(betas_second_stage_estimated_kyri_1(1,:))
beta2_second_stage_estimated_kyri_1 = mean(betas_second_stage_estimated_kyri_1(2,:))-2.5
beta2_std_second_stage_estimated_kyri_1 = std(betas_second_stage_estimated_kyri_1(2,:))

betas_second_stage_estimated_kyri_2 = betas_second_stage_estimated_kyri(:,3);
betas_second_stage_estimated_kyri_2 = reshape(betas_second_stage_estimated_kyri_2, 2,500);
beta1_second_stage_estimated_kyri_2 = mean(betas_second_stage_estimated_kyri_2(1,:))-1
beta1_std_second_stage_estimated_kyri_2 = std(betas_second_stage_estimated_kyri_2(1,:))
beta2_second_stage_estimated_kyri_2 = mean(betas_second_stage_estimated_kyri_2(2,:))-2.5
beta2_std_second_stage_estimated_kyri_2 = std(betas_second_stage_estimated_kyri_2(2,:))

betas_second_stage_estimated_kyri_3 = betas_second_stage_estimated_kyri(:,4);
betas_second_stage_estimated_kyri_3 = reshape(betas_second_stage_estimated_kyri_3, 2,500);
beta1_second_stage_estimated_kyri_3 = mean(betas_second_stage_estimated_kyri_3(1,:))-1
beta1_std_second_stage_estimated_kyri_3 = std(betas_second_stage_estimated_kyri_3(1,:))
beta2_second_stage_estimated_kyri_3 = mean(betas_second_stage_estimated_kyri_3(2,:))-2.5
beta2_std_second_stage_estimated_kyri_3 = std(betas_second_stage_estimated_kyri_3(2,:))

betas_second_stage_estimated_kyri_5 = betas_second_stage_estimated_kyri(:,5);
betas_second_stage_estimated_kyri_5 = reshape(betas_second_stage_estimated_kyri_5, 2,500);
save('/Users/gabrielaszini/Documents/Master Thesis/Simulations/Figures/beta_5_kyriazidou_Design5.mat','betas_second_stage_estimated_kyri_5')
beta1_second_stage_estimated_kyri_5 = mean(betas_second_stage_estimated_kyri_5(1,:))-1
beta1_std_second_stage_estimated_kyri_5 = std(betas_second_stage_estimated_kyri_5(1,:))
beta2_second_stage_estimated_kyri_5 = mean(betas_second_stage_estimated_kyri_5(2,:))-2.5
beta2_std_second_stage_estimated_kyri_5 = std(betas_second_stage_estimated_kyri_5(2,:))

betas_second_stage_estimated_kyri_10 = betas_second_stage_estimated_kyri(:,6);
betas_second_stage_estimated_kyri_10 = reshape(betas_second_stage_estimated_kyri_10, 2,500);
save('/Users/gabrielaszini/Documents/Master Thesis/Simulations/Figures/beta_10_kyriazidou_Design5.mat','betas_second_stage_estimated_kyri_10')
beta1_second_stage_estimated_kyri_10 = mean(betas_second_stage_estimated_kyri_10(1,:))-1
beta1_std_second_stage_estimated_kyri_10 = std(betas_second_stage_estimated_kyri_10(1,:))
beta2_second_stage_estimated_kyri_10 = mean(betas_second_stage_estimated_kyri_10(2,:))-2.5
beta2_std_second_stage_estimated_kyri_10 = std(betas_second_stage_estimated_kyri_10(2,:))

% with estimated first stage and corrected
betas_second_stage_estimated_corrected_kyri = cell2mat(OUTPUT(:,10));
betas_second_stage_estimated_corrected_kyri_05 = betas_second_stage_estimated_corrected_kyri(:,1);
betas_second_stage_estimated_corrected_kyri_05 = reshape(betas_second_stage_estimated_corrected_kyri_05, 2,500);
beta1_second_stage_estimated_corrected_kyri_05 = mean(betas_second_stage_estimated_corrected_kyri_05(1,:))-1
beta1_std_second_stage_estimated_corrected_kyri_05 = std(betas_second_stage_estimated_corrected_kyri_05(1,:))
beta2_second_stage_estimated_corrected_kyri_05 = mean(betas_second_stage_estimated_corrected_kyri_05(2,:))-2.5
beta2_std_second_stage_estimated_corrected_kyri_05 = std(betas_second_stage_estimated_corrected_kyri_05(2,:))

betas_second_stage_estimated_corrected_kyri_1 = betas_second_stage_estimated_corrected_kyri(:,2);
betas_second_stage_estimated_corrected_kyri_1 = reshape(betas_second_stage_estimated_corrected_kyri_1, 2,500);
beta1_second_stage_estimated_corrected_kyri_1 = mean(betas_second_stage_estimated_corrected_kyri_1(1,:))-1
beta1_std_second_stage_estimated_corrected_kyri_1 = std(betas_second_stage_estimated_corrected_kyri_1(1,:))
beta2_second_stage_estimated_corrected_kyri_1 = mean(betas_second_stage_estimated_corrected_kyri_1(2,:))-2.5
beta2_std_second_stage_estimated_corrected_kyri_1 = std(betas_second_stage_estimated_corrected_kyri_1(2,:))

betas_second_stage_estimated_corrected_kyri_2 = betas_second_stage_estimated_corrected_kyri(:,3);
betas_second_stage_estimated_corrected_kyri_2 = reshape(betas_second_stage_estimated_corrected_kyri_2, 2,500);
beta1_second_stage_estimated_corrected_kyri_2 = mean(betas_second_stage_estimated_corrected_kyri_2(1,:))-1
beta1_std_second_stage_estimated_corrected_kyri_2 = std(betas_second_stage_estimated_corrected_kyri_2(1,:))
beta2_second_stage_estimated_corrected_kyri_2 = mean(betas_second_stage_estimated_corrected_kyri_2(2,:))-2.5
beta2_std_second_stage_estimated_corrected_kyri_2 = std(betas_second_stage_estimated_corrected_kyri_2(2,:))

betas_second_stage_estimated_corrected_kyri_3 = betas_second_stage_estimated_corrected_kyri(:,4);
betas_second_stage_estimated_corrected_kyri_3 = reshape(betas_second_stage_estimated_corrected_kyri_3, 2,500);
beta1_second_stage_estimated_corrected_kyri_3 = mean(betas_second_stage_estimated_corrected_kyri_3(1,:))-1
beta1_std_second_stage_estimated_corrected_kyri_3 = std(betas_second_stage_estimated_corrected_kyri_3(1,:))
beta2_second_stage_estimated_corrected_kyri_3 = mean(betas_second_stage_estimated_corrected_kyri_3(2,:))-2.5
beta2_std_second_stage_estimated_corrected_kyri_3 = std(betas_second_stage_estimated_corrected_kyri_3(2,:))

betas_second_stage_estimated_corrected_kyri_5 = betas_second_stage_estimated_corrected_kyri(:,5);
betas_second_stage_estimated_corrected_kyri_5 = reshape(betas_second_stage_estimated_corrected_kyri_5, 2,500);
beta1_second_stage_estimated_corrected_kyri_5 = mean(betas_second_stage_estimated_corrected_kyri_5(1,:))-1
beta1_std_second_stage_estimated_corrected_kyri_5 = std(betas_second_stage_estimated_corrected_kyri_5(1,:))
beta2_second_stage_estimated_corrected_kyri_5 = mean(betas_second_stage_estimated_corrected_kyri_5(2,:))-2.5
beta2_std_second_stage_estimated_corrected_kyri_5 = std(betas_second_stage_estimated_corrected_kyri_5(2,:))

betas_second_stage_estimated_corrected_kyri_10 = betas_second_stage_estimated_corrected_kyri(:,6);
betas_second_stage_estimated_corrected_kyri_10 = reshape(betas_second_stage_estimated_corrected_kyri_10, 2,500);
beta1_second_stage_estimated_corrected_kyri_10 = mean(betas_second_stage_estimated_corrected_kyri_10(1,:))-1
beta1_std_second_stage_estimated_corrected_kyri_10 = std(betas_second_stage_estimated_corrected_kyri_10(1,:))
beta2_second_stage_estimated_corrected_kyri_10 = mean(betas_second_stage_estimated_corrected_kyri_10(2,:))-2.5
beta2_std_second_stage_estimated_corrected_kyri_10 = std(betas_second_stage_estimated_corrected_kyri_10(2,:))


