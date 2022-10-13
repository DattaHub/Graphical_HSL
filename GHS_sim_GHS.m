%%%% estimate precision matrix by Graphical Horseshoe
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
    [GHS_omega_save,GHS_omega_vector_save,~,~] = GHS(S,n,500,1000); %burnin = 100 and nmc = 5000 in paper 
	GHS_time(i) = cputime-t;	
    GHS_est(:,:,i) = mean(GHS_omega_save,3);    % estimate of precision matrix
	% credible interval
	a = 50;
	for j = 1:(p*(p-1)/2)
		GHS_omega_vector_lb(i,j) = prctile(GHS_omega_vector_save(j,:),(100-a)/2);
		GHS_omega_vector_ub(i,j) = prctile(GHS_omega_vector_save(j,:),100-(100-a)/2);
		GHS_omega_zero(i,j) = GHS_omega_vector_lb(i,j)<0 & GHS_omega_vector_ub(i,j)>0;
    end
end

for i = 1:m
    GHS_omega_est = GHS_est(:,:,i);
    GHS_sigma_est = inv(GHS_omega_est);
    % Stein's loss of Sigma_inv; Frobenius norm of sigma_inv-Sigma_inv
    GHS_Steinsloss(i) = log(det(GHS_sigma_est*sigma_inv))+trace(GHS_omega_est*sigma)-p;
    GHS_Fnorm(i) = norm(GHS_omega_est-sigma_inv,'fro');
end

for i = 1:m
	TP = sum((GHS_omega_zero(i,:)==0).*(omega_elements~=0));
	FP = sum((GHS_omega_zero(i,:)==0).*(omega_elements==0));
	TN = sum((GHS_omega_zero(i,:)~=0).*(omega_elements==0));
	FN = sum((GHS_omega_zero(i,:)~=0).*(omega_elements~=0));
    
    MCC_matrix(i) = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    
	sen_GHS(i) = TP/sum(omega_elements~=0);    % true positive rate
    spe_GHS(i) = TN/sum(omega_elements==0);
	fpr_GHS(i) = 1-spe_GHS(i); % false positive rate 
end

means = [mean(GHS_Steinsloss),mean(GHS_Fnorm),mean(sen_GHS),mean(fpr_GHS), mean(MCC_matrix), mean(GHS_time)];
fprintf('GHS means: loss, Fnorm, TPR, FPR, MCC, time %f, %f, %f, %f, %f, %f',means)

sds = [std(GHS_Steinsloss),std(GHS_Fnorm),std(sen_GHS),std(fpr_GHS), std(MCC_matrix), std(GHS_time)];
fprintf('GHS sds: loss, Fnorm, TPR, FPR %f, %f, %f, %f, %f, %f',sds)