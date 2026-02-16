
function [dF,n_opt,idx_opt,p_perm,qFDR,SES,...
    Num_nodes,Npd_mean,Npd_sigma,tau,KE,freq] = ...
    partition_recur(idx_cell,J,Cp,pC,flag_rndGraph,Nperm,n_large,n_min)

for i = 1:size(idx_cell,1)
    idx0 = idx_cell(i,:);
    k = max(idx0); % k = #partitions = 2, except first scale (k=1)

    if isequal(idx0,ones(size(idx0))) || (k==2) % skip indivisible partition (terminated branch)
        
        for j = 1:k % for each partition
            idx_range = find(idx0 == j);
            n = length(idx_range); % n = network size
            if length(idx_range)>= n_min % min size of graph (#nodes) to partition
                [dF(i,j),n_opt(i,j),idx_opt(i,j,:),p_perm(i,j),SES(i,j),...
                    Num_nodes(i,j),Npd_mean(i,j),Npd_sigma(i,j),...
                    tau{i,j},KE{i,j},freq{i,j}] = ...
                    partition(idx_range,J,Cp,pC,flag_rndGraph,Nperm,n_large);
            else
                dF(i,j) = 0;
                n_opt(i,j) = 1;
                idx_opt(i,j,1:size(idx_cell,2)) = 0;
                idx_opt(i,j,idx_range) = 1; % (1's for indivis partition)
                NCut(i,j) = 0;
                p_perm(i,j) = 1; % indivisible
                SES(i,j) = -8; % dummy flag
                Num_nodes(i,j) = length(idx_range);
                Npd_mean(i,j) = 0;
                Npd_sigma(i,j) = 0;
                tau{i,j} = 0;
                KE{i,j} = 0;
                freq{i,j}=0;
            end
        end
    end
end

%% Remove fully terminated branches, 
% and perform family-wise FDR correction at each scale
qFDR = p_perm; % initial
if exist('n_opt','var')
    % find and remove indivisible partitions from previous scale (ie terminated branches) 
    % that were *not* tested at this scale
    null = any(n_opt==0,2); 
    dF(null,:) = []; % remove
    n_opt(null, :) = [];
    idx_opt(null,:,:) = [];
    p_perm(null,:) = [];
    SES(null,:) = [];
    Num_nodes(null,:) = [];
    Npd_mean(null,:) = [];
    Npd_sigma(null,:) = [];
    tau(null,:) = [];
    KE(null,:) = [];
    freq(null,:) = [];

    % FDR control over all families of hypotheses at each scale;
    % Refs: (Yekutieli,2006, pages 7-8) and (Benjamini & Bogomolov,2014, Eq. 3)
    % Controlling FDR <= q for each family, guarantees avg FDR <= q over all families 
    qFDR = p_perm;
    for i = 1:size(p_perm,1) % for each family of hypotheses
        p = p_perm(i,:); % each row = a family
        if ~isempty(p)
            qFDR (i,:) = mafdr(p,'BHFDR',true);
            n_opt(i,:) = (qFDR(i,:) <0.05) + 1; % set n_opt = 1 if qFDR >= 0.05
            idx_opt(i,(n_opt(i,:)==1),:) = ... 
                logical(idx_opt(i,(n_opt(i,:)==1),:)); % (1s for indivisible partition based on qFDR)
        end
    end
    % since indivisible partitions were already included in scale analysis of prev scale:
    tau(n_opt==1) = {0} ; 
    KE(n_opt==1)  = {0} ; 
    freq(n_opt==1)= {0} ; 
end
