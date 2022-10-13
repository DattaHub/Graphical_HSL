
n_EMs = 50; % 50 denotes the 50 starting points with which we start the 
% MAP estimation with horseshoe-like penalty

% starting point generation
p = 100;
%p = 200;
Omega_saves = zeros(p,p,n_EMs);
rng(123456789);

for i = 1:n_EMs %choose parfor instead of "for" loop if running on a multi core CPU

start_point = eye(p);

    for row = 2:p
        d = 0;
            while d~=p
                row_seq = row:1:p;
                col_seq = 1:1:(p-row+1);

                rand_noise = -0.05 + rand(1, length(row_seq))*2*0.05;
                
                lin_idcs = sub2ind(size(start_point), row_seq, col_seq);
                start_point(lin_idcs) = rand_noise;

                lin_idcs = sub2ind(size(start_point), col_seq, row_seq);
                start_point(lin_idcs) = rand_noise;

                d = eig(start_point);
                d = sum(d>0);

            end 
    end 
    
fprintf("Finished %d th start point generation out of %d start points \n", i, n_EMs);
Omega_saves(:,:,i) = start_point;
end 

total_sets = 1;
%%%%% total_sets = 50 in paper

n = 120;

precision_str = 'randomn';
%precision_str = 'hubs';
%precision_str = 'cliquespos75';
%precision_str = 'cliquesneg45';

for file_iter = 1:total_sets
    
    fprintf('Data set number - %d is being processed \n ',file_iter);
    FileName=['GHS_sim_p',num2str(p),precision_str,num2str(n),'_data',num2str(file_iter),'.csv'];
    xx = csvread(FileName);
    [n,q] = size(xx);
    S = xx'*xx;

    [Omega_est, final_nu_matrix, total_iterations, each_time_taken] =...
        Multi_start_point_Fixed_b_EM_HS_like(Omega_saves, S, n, q, n_EMs);
    
    FileName=['Workspace_of',num2str(file_iter),'th_data_set_in',num2str(n_EMs),...
        'EMs','.mat'];
    
    save(FileName, 'Omega_est', 'final_nu_matrix', 'total_iterations','each_time_taken');
    fprintf('Data set number - %d is finished \n ',file_iter);
    
end 

% run the next few lines to get Stein's loss, Forbenious Norm %%%%
% TPR, FPR, MCC and time taken %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Omega_each_DS = zeros(q,q,total_sets);

for i = 1:total_sets
    
    file_iter = i;
    FileName=['Workspace_of',num2str(file_iter),'th_data_set_in',num2str(n_EMs),...
        'EMs','.mat'];
   
    matObj = matfile(FileName);
    
    Omega_each_DS(:,:,i) = mean(matObj.Omega_est,3);
end 

True_Omega = readmatrix(['GHS_sim_p',num2str(p),precision_str,'_sigmainv.csv']);

omega_elements = True_Omega(tril(true(size(True_Omega)),-1))'; 

tp_tn_fp_fn_matrix = zeros(total_sets, 4);

for i = 1:total_sets
    
    temp_array = Omega_each_DS(:,:,i);
    omega_elements_current = temp_array(tril(true(size(temp_array)),-1))'; 
    % true Positive    
    tp_tn_fp_fn_matrix(i,1) = sum(omega_elements~=0 & omega_elements_current~=0);
    % true negative    
    tp_tn_fp_fn_matrix(i,2) = sum(omega_elements==0 & omega_elements_current==0);
    % false positive     
    tp_tn_fp_fn_matrix(i,3) = sum(omega_elements==0 & omega_elements_current~=0);
    % false negative    
    tp_tn_fp_fn_matrix(i,4) = sum(omega_elements~=0 & omega_elements_current==0);
end 

tpr_fpr_matrix = zeros(total_sets,2);

for i = 1: total_sets
   
        tpr_fpr_matrix(i,1) = tp_tn_fp_fn_matrix(i,1)/sum(omega_elements~=0);
        tpr_fpr_matrix(i,2) = 1- (tp_tn_fp_fn_matrix (i,2)/sum(omega_elements==0));
end 

MCC_matrix = zeros(1,total_sets);

for i = 1:total_sets
    TP = tp_tn_fp_fn_matrix(i,1);
    TN = tp_tn_fp_fn_matrix(i,2);
    FP = tp_tn_fp_fn_matrix(i,3);
    FN = tp_tn_fp_fn_matrix(i,4);
    
    MCC_matrix(i) = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
end 

diff_Frobenious_norm = zeros(1, total_sets);

for i = 1:total_sets
    
   diff_Frobenious_norm(i) = norm(Omega_each_DS(:,:,i) - True_Omega, 'fro');
end 

stein_loss = zeros(1, total_sets);

for i = 1:total_sets
    
    stein_loss(i) = -log(det(True_Omega\Omega_each_DS(:,:,i)))+...
        trace(True_Omega\Omega_each_DS(:,:,i))-q;
    
end 

time_taken = zeros(1, total_sets);

for i = 1:total_sets
    
    file_iter = i;
    FileName=['Workspace_of',num2str(file_iter),'th_data_set_in',num2str(n_EMs),...
        'EMs','.mat'];
   
    matObj = matfile(FileName);
    
    time_taken(i)= mean(matObj.each_time_taken);
    
end


fprintf("Mean Stein's Loss is,  %f, std of stein loss is  %f \n" ,round(mean(stein_loss),3), ...
    round(std(stein_loss ),3));

fprintf("Mean Frobenious norm is,  %f, std of From Norm is  %f \n" ,round(mean(diff_Frobenious_norm),3), ...
    round(std(diff_Frobenious_norm),3));

fprintf("Mean TPR, is,  %f, std of TPR is %f \n " ,round(mean(tpr_fpr_matrix(:,1)),3), ...
    round(std(tpr_fpr_matrix(:,1) ),3));

fprintf("Mean FPR is,  %f, std of FPR is  %f\n " ,round(mean(tpr_fpr_matrix(:,2),1),3), ...
    round(std(tpr_fpr_matrix(:,2) ),3));

fprintf("Mean MCC is,  %f, std of MCC  %f \n" ,round(mean(MCC_matrix),3), ...
    round(std(MCC_matrix),3));

fprintf("Mean Time taken is,  %f, std of Time taken is  %f \n" ,round(mean(time_taken),3), ...
    round(std(time_taken),3));