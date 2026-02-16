
function [dF,n_opt,idx_opt,p_perm,SES,n,Npd_mean,Npd_sigma,tau,KE,freq] = ...
    partition(idx_range,J0,Cp0,pC0,flag_rndGraph,Nperm,n_large)

n0 = size(J0,1);
if isempty(idx_range)==1, idx_range=1:n0;end
J = J0(idx_range,idx_range);
n = size(J,1); % #nodes/regions
Cp = Cp0(idx_range,idx_range); % posterior variance
pC = pC0(idx_range,idx_range); % prior variance
k_partitions = 2; % #partitions
p_perm = 1; % initial value
Npd_mean = 0; % avg minCut
Npd_sigma = 0; % std of minCut
tau = 0; % dissipative time constant
KE = 0; % kinetic energy
freq = 0; % instrisic freq
SES = -10;  % dummy initial value (when not permuting)
partition_evid = 0; % in general will run permutation
if exist('Nperm','var') && numel(Nperm)==2
    if n>n_large, Nperm = Nperm(1);
    else, Nperm = Nperm(2);
    end
end
if Nperm==0, partition_evid = 1; p_perm = 0; end % assume divisible partition

%% Spectral clustering (SC) on transformed graph

% Cp = reshape(diag(DCM.Cp),n,n); So, Cp is n*n (not n^2 * n^2)
Pp = 1./Cp; 
J_tf = (Pp).*(J.^2) + log(Cp) - log(pC);
J_tf (isnan(J_tf) | J_tf==Inf | J_tf ==-Inf) = 0;
J_tf_sym = (J_tf + J_tf')/2;
lap_type = 'normalized';

if flag_rndGraph == 1
    [J_tf_sym_rnd] = randmio_dir_connected(J_tf_sym,1); % BCT toolbox
    J_tf_sym = (J_tf_sym_rnd + J_tf_sym_rnd')/2;
end

% Apply SC to J_tf_sym & compute dF
G = [];
param_CSC.lap_type = lap_type;  % type of Laplacian used; 'normalized' or 'combinatorial'
G.N = n;                        % number of nodes
G.W = J_tf_sym;                 % adjacanecy
G.k = k_partitions;             % # partitions for SC
[IDX_J] = SC_filloutlier(G,param_CSC.lap_type,1);

if flag_rndGraph == 0 % --> Compute dF using nBMR & J
    % construct block-diagonal J_pt from J
    J_pt = zeros(n);
    for i = 1:k_partitions % for each partition
        pt = zeros(n,n);
        range = (IDX_J == i);
        pt(range,range) = J(range,range);
        J_pt = pt + J_pt;
    end
    % Compute dF using nBMR & J
    [dF,sE,sC] = spm_log_evidence_naive(J,J_pt,Cp,pC);

else % if flag_rndGraph == 1 --> Compute dF = -Cut
    % construct block-diagonal J_pt from J_tf_sym
        J_pt = zeros(n);
    for i = 1:k_partitions % for each partition
        pt = zeros(n,n);
        range = (IDX_J == i);
        pt(range,range) = J_tf_sym(range,range); 
        J_pt = pt + J_pt;
    end
    % Compute dF = -Cut
    dF = -0.5*sum(J_tf_sym - J_pt,'all'); % dF = -Cut= -0.5* sum_of_off_diag
end

%% Graph randomization (for partition significance testing)

if partition_evid == 0 && flag_rndGraph == 0
    flag_rndGraph = 1;
    for perm_iter = 1:Nperm
    [dFmax_rnd(perm_iter)] = partition(idx_range,J0,Cp0,pC0,flag_rndGraph,[],n_large);
    perm_iter
    end
    flag_rndGraph = 0;
    p_perm = (sum(dFmax_rnd >= dF)+1)/(Nperm+1); % conservative p-value
    partition_evid = (p_perm < 0.05);
    SES = (dF - mean(dFmax_rnd))/std(dFmax_rnd); % standardized effect size

    % Schreiber & Martin,1999: self-averaging property of minCut
    % --> std(minCut)/mean(minCut) --> 0 as N --> Inf
    minCut_rnd = -dFmax_rnd';
    Npd = fitdist(minCut_rnd,'Normal'); 
    Npd_sigma = std(Npd);
    Npd_mean = mean(Npd);
end

%% Compute dynamical measures (Tau, KE, freq)
if partition_evid ==1 && flag_rndGraph == 0
    run compute_dyn_meas.m
end

%% Prepare output (n_opt and idx_opt)
if partition_evid ==1 && flag_rndGraph == 0
    n_opt = 2 ; % divisible partition
    idx_opt_local = IDX_J;
    idx_opt = zeros(1,n0);idx_opt(idx_range)=idx_opt_local;

else % if partition_evid == 0 || flag_rndGraph == 1
    n_opt = 1; % indivisible partition
    idx_opt = zeros(1,n0); idx_opt(idx_range)= ones(size(idx_range));
end



