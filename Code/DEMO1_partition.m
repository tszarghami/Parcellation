
%% Multiscale Parcellation of a DCM 
% from attention-to-visual-motion fMRI datatset
% Tahereh Zarghami, 2025/6/15

%% Load or compute effective connetivity matrix (J)

clear;clc;close all

% Get current script's directory and go to root directory
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
cd(rootDir)
addpath(genpath(rootDir));

% Load the precomputed effective connectivity matrix (J),
% which was computed by applying linear DCM to the attention dataset
load('J_1024.mat')
% Alternatively, to compute J by yourself, follow compute_J.m

% Results directory
mkdir('Results-DEMO1')
cd('Results-DEMO1')

%% Recursive bipartitioning

n = size(J,1);
IDX_cell{1} = ones(1,n);
flag_continue = 1;
scale = 1;
flag_rndGraph = 0;
n_large = 100; % define large graphs > 100 
n_min = 8; % min graph size to partition
Nperm = [0,100]; % # permutations (i.e. randomized graphs)
% first element = #perm for large graphs
% since large graphs of this dataset were divisible,
% we bypass permutations for sake of demo time

while flag_continue 

    [dF{scale},N_opt{scale},IDX_raw{scale},P_perm{scale}, ...
        qFDR{scale},SES{scale},Num_nodes{scale},...
        Npd_Mean{scale},Npd_Sigma{scale},...
        Tau{scale},KE{scale},Freq{scale}] = ...
        partition_recur(IDX_cell{scale},J,Cp,pC,flag_rndGraph,Nperm,n_large,n_min);

    idx = IDX_raw{scale};
    idx_stacked = [];
    for  k = 1:size(idx,1)
        idx_sqz = squeeze(idx(k,:,:));
        if iscolumn(idx_sqz), idx_sqz = idx_sqz';end
        idx_stacked = [idx_stacked; idx_sqz];
    end
    flag_continue = any(any(N_opt{scale}>1));

    scale = scale+1;
    IDX_cell{scale} = idx_stacked;
    
end

%% Spatial maps of partitions at different scales (as MIPs)

scales_to_plot = [1:5]; % specify the scales to visualize
max_partitions = 100; % specifiy max #partitions from each scale to visualize
run plot_anatomy.m 

%% Scale invariance (log-log plots)

run plot_loglog.m

