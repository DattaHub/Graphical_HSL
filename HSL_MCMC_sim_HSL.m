% estimate precision matrix by Graphical Horseshoe like
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
    S = xx'*xx;
	t = cputime;
    [HSL_MC_omega_save,HSL_MC_omega_vector_save,~,~] = HSL_MCMC(S,n,500,1000); %burnin = 1000 and nmc = 5000 in paper 
	HSL_MC_time(i) = cputime-t;	
    HSL_MC_est(:,:,i) = mean(HSL_MC_omega_save,3);    % estimate of precision matrix
	% credible interval
	a = 50;
	for j = 1:(p*(p-1)/2)
		HSL_MC_omega_vector_lb(i,j) = prctile(HSL_MC_omega_vector_save(j,:),(100-a)/2);
		HSL_MC_omega_vector_ub(i,j) = prctile(HSL_MC_omega_vector_save(j,:),100-(100-a)/2);
		HSL_MC_omega_zero(i,j) = HSL_MC_omega_vector_lb(i,j)<0 & HSL_MC_omega_vector_ub(i,j)>0;
    end
end 

for i = 1:m
    HSL_MC_omega_est = HSL_MC_est(:,:,i);
    HSL_MC_sigma_est = inv(HSL_MC_omega_est);
    % Stein's loss of Sigma_inv; Frobenius norm of sigma_inv-Sigma_inv
    HSL_MC_Steinsloss(i) = log(det(HSL_MC_sigma_est*sigma_inv))+trace(HSL_MC_omega_est*sigma)-p;
    HSL_MC_Fnorm(i) = norm(HSL_MC_omega_est-sigma_inv,'fro');
end

for i = 1:m
	TP = sum((HSL_MC_omega_zero(i,:)==0).*(omega_elements~=0));
	FP = sum((HSL_MC_omega_zero(i,:)==0).*(omega_elements==0));
	TN = sum((HSL_MC_omega_zero(i,:)~=0).*(omega_elements==0));
	FN = sum((HSL_MC_omega_zero(i,:)~=0).*(omega_elements~=0));
    
    MCC_matrix(i) = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    
	sen_HSL_MC(i) = TP/sum(omega_elements~=0);    % true positive rate
    spe_HSL_MC(i) = TN/sum(omega_elements==0);
	fpr_HSL_MC(i) = 1-spe_HSL_MC(i); % false positive rate 
end

means = [mean(HSL_MC_Steinsloss),mean(HSL_MC_Fnorm),mean(sen_HSL_MC),mean(fpr_HSL_MC), mean(MCC_matrix), mean(HSL_MC_time)];
fprintf('HSL_MC means: loss, Fnorm, TPR, FPR, MCC, time %f, %f, %f, %f, %f, %f',means)

sds = [std(HSL_MC_Steinsloss),std(HSL_MC_Fnorm),std(sen_HSL_MC),std(fpr_HSL_MC), std(MCC_matrix), std(HSL_MC_time)];
fprintf('HSL_MC sds: loss, Fnorm, TPR, FPR %f, %f, %f, %f, %f, %f',sds)