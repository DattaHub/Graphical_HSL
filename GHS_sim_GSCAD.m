% Estimate precision matrix by GSCAD, using Hao Wang (2012)'s code

precision_str = 'randomn';
%precision_str = 'hubs';
%precision_str = 'cliquespos75';
%precision_str = 'cliquesneg45';

p = 100;
% p= 200;
n = 120;

sigma_inv = csvread(['GHS_sim_p',num2str(p),precision_str,'_sigmainv.csv']);
sigma = inv(sigma_inv);
omega_elements = sigma_inv(tril(true(size(sigma_inv)),-1))';    % elements in the upper triangular sigma_inv

m = 1;
% m = 50 in paper 
for i = 1:m
    fprintf('Data set number - %d is being processed \n ',i);
    FileName=['GHS_sim_p',num2str(p),precision_str,num2str(n),'_data',num2str(i),'.csv'];
    xx = csvread(FileName);
    [n,p] = size(xx);

    Y = xx';
    rhopath = exp(-3:0.1:0);    % path of parameter rho; tuning parameter a is fixed at 3.7
	t = cputime;
    [GSCAD_sigma_est,GSCAD_omega_est,GSCAD_M0,GSCAD_rhomax] = glasso_SCAD_cv(Y,rhopath,5);
	GSCAD_time(i) = cputime-t;

    GSCAD_est(:,:,i) = GSCAD_omega_est;    % estimate of precision matrix
    GSCAD_cov_est(:,:,i) = GSCAD_sigma_est;   % estimate of covariance matrix
    GSCAD_rho(i) = GSCAD_rhomax;    % tuning parameter rho used in estimate
end

for i = 1:m
    GSCAD_omega_est = GSCAD_est(:,:,i);
    GSCAD_sigma_est = GSCAD_cov_est(:,:,i);   
    % Stein's loss of Sigma_inv; Frobenius norm of sigma_inv-Sigma_inv
    GSCAD_Steinsloss(i) = log(det(GSCAD_sigma_est*sigma_inv))+trace(GSCAD_omega_est*sigma)-p;
    GSCAD_Fnorm(i) = norm(GSCAD_omega_est-sigma_inv,'fro');

    GSCAD_omega_elements = GSCAD_omega_est(tril(true(size(sigma_inv)),-1))';
	TP = sum((abs(GSCAD_omega_elements)~=0).*(omega_elements~=0));
	TN = sum((abs(GSCAD_omega_elements)==0).*(omega_elements==0));
	FP = sum((abs(GSCAD_omega_elements)~=0).*(omega_elements==0));
    FN = sum((abs(GSCAD_omega_elements)==0).*(omega_elements~=0));
    MCC_matrix(i) = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    GSCAD_sen(i) = TP/sum(omega_elements~=0);    % true positives rate
	GSCAD_spe(i) = TN/sum(omega_elements==0);
    GSCAD_fpr(i) = 1-GSCAD_spe(i);    % false positive rate
end

means = [mean(GSCAD_Steinsloss),mean(GSCAD_Fnorm),mean(sen_GSCAD),mean(fpr_GSCAD), mean(MCC_matrix), mean(GSCAD_time)];
fprintf('GSCAD means: loss, Fnorm, TPR, FPR, MCC, time %f, %f, %f, %f, %f, %f',means)

sds = [std(GSCAD_Steinsloss),std(GSCAD_Fnorm),std(sen_GSCAD),std(fpr_GSCAD), std(MCC_matrix), std(GSCAD_time)];
fprintf('GSCAD sds: loss, Fnorm, TPR, FPR %f, %f, %f, %f, %f, %f',sds)