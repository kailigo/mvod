clear;
close all;

addpath data

outlier_ratios = [0.05, 0.05, 0.05];
nview = 2;

lambda1 = 0.01;
alpha = 1;
beta = 1;

n = 20;

parfor i =  1:n    
    [view_data, classLabel, outLabel] =  gen_outliers_zoo( outlier_ratios, nview );    
    AUC_val = LDSR(view_data, classLabel, outLabel, lambda1, alpha, beta);     
    aucs(i, :) =  AUC_val;
end
mean_aucs = mean(aucs);
[max_aucs, idx] = max(mean_aucs);
auc = aucs(:, idx);

fprintf('Mean AUC is: %f\n', mean(auc));
fprintf('AUC std is: %f\n', std(auc));
