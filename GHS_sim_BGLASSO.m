% Bayesian graphical lasso estimate of precision matrix
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
% Bayesian graphical lasso estimate, using Hao Wang (2012)'s code

for i=1:m
    fprintf('Data set number - %d is being processed \n ',i);
    FileName=['GHS_sim_p',num2str(p),precision_str,num2str(n),'_data',num2str(i),'.csv'];
    xx = csvread(FileName);
    [n,p] = size(xx);
    S = xx'*xx;
	t = cputime;
	% use identity matrix as starting values
	% use Gamma(1,0.01) as hyperprior of \lambda, as in Hao Wang (2012)
	[BGLASSO_sigma_save,BGLASSO_omega_save,~,BGLASSO_omega_vector_save] = ...
		BayesGLasso_Columnwise(S,n,eye(p),eye(p),1,0.01,500,1000);
	BGLASSO_time(i) = cputime-t;

    BGLASSO_est(:,:,i) = mean(BGLASSO_omega_save,3);    % estimate of precision matrix
    BGLASSO_cov_est(:,:,i) = mean(BGLASSO_sigma_save,3);   % estimate of covariance matrix
	a = 50;
	for j = 1:(p*(p-1)/2)
		BGLASSO_omega_vector_lb(i,j) = prctile(BGLASSO_omega_vector_save(j,:),(100-a)/2);
		BGLASSO_omega_vector_ub(i,j) = prctile(BGLASSO_omega_vector_save(j,:),100-(100-a)/2);
		BGLASSO_omega_zero(i,j) = BGLASSO_omega_vector_lb(i,j)<0 & BGLASSO_omega_vector_ub(i,j)>0;
    end
end

for i=1:m
    BGLASSO_omega_est = BGLASSO_est(:,:,i);
    BGLASSO_sigma_est = BGLASSO_cov_est(:,:,i);
    % Stein's loss of Sigma_inv; Frobenius norm of sigma_inv-Sigma_inv
    BGLASSO_Steinsloss(i) = log(det(BGLASSO_sigma_est*sigma_inv))+trace(BGLASSO_omega_est*sigma)-p;
    BGLASSO_Fnorm(i) = norm(BGLASSO_omega_est-sigma_inv,'fro');
end

for i = 1:m
	TP = sum((BGLASSO_omega_zero(i,:)==0).*(omega_elements~=0));
	FP = sum((BGLASSO_omega_zero(i,:)==0).*(omega_elements==0));
	TN = sum((BGLASSO_omega_zero(i,:)~=0).*(omega_elements==0));
	FN = sum((BGLASSO_omega_zero(i,:)~=0).*(omega_elements~=0));
    MCC_matrix(i) = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
	sen_BGLASSO(i) = TP/sum(omega_elements~=0);    % true positive rate
	spe_BGLASSO(i) = TN/sum(omega_elements==0);
	fpr_BGLASSO(i) = 1-spe_BGLASSO(i); % False positive rate 
end

means = [mean(BGLASSO_Steinsloss),mean(BGLASSO_Fnorm),mean(sen_BGLASSO),mean(fpr_BGLASSO), mean(MCC_matrix), mean(BGLASSO_time)];
fprintf('BGLASSO means: loss, Fnorm, TPR, FPR, MCC, time %f, %f, %f, %f, %f, %f',means)

sds = [std(BGLASSO_Steinsloss),std(BGLASSO_Fnorm),std(sen_BGLASSO),std(fpr_BGLASSO), std(MCC_matrix), std(BGLASSO_time)];
fprintf('BGLASSO sds: loss, Fnorm, TPR, FPR %f, %f, %f, %f, %f, %f',sds)