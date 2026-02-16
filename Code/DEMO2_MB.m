
%% DEMO: Identify MB based on effective connectivity
% and assess its validity based on functional connectivity
% Tahereh Zarghami, 2025/6/15

%% Compute H = Inverse Cov from empirical J
% clear;clc;close all

% Get current script's directory and go to root directory
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
cd(rootDir)
addpath(genpath(rootDir));

load('J_71.mat')
mkdir('Results-DEMO2')
cd('Results-DEMO2')

n = size(J,1); % #nodes
Gamma = eye(n); % since P{1}.h = 1 approx in spm_dcm_J.m
H = spm_ness_mod (J,Gamma*eye(n));
Nsample = 360; % #time points in attention dataset
n = size(J,1); % Choose J s.t. n < Nsample

%% Compute partial correlation (pcorr) from H

% Compute partial correlation from H
Hdiag = diag(H);
pcorr = zeros(n,n);
for i = 1:size(H,1)
    for j = i:size(H,2)
    pcorr (i,j) = -H(i,j)/sqrt(Hdiag(i)*Hdiag(j));
    end
end
pcorr = pcorr + pcorr' - diag(diag(pcorr));

%% Bayesian hypothesis testing on partial correlation

% Ref: Kucharský, Š., Wagenmakers, E. J., van den Bergh, D., & Ly, A. (2023). 
% Analytic Posterior Distribution and Bayes Factor for Pearson Partial Correlations.
p = n-2;% #variables in condition set
alpha = 1/2; % alpha = hyperparameter of stretched beta prior
alpha_tilda = alpha + (Nsample-p-1)/2;
a = (Nsample-p-1)/2;
b = (Nsample-p-1)/2;
c = alpha_tilda+ 1/2;
BF10 = zeros(n);
for i = 1:size(pcorr,1)
    for j = i:size(pcorr,2)
    pc = pcorr(i,j);
    BF10 (i,j) = beta(1/2,alpha_tilda)/beta(1/2,alpha)*hypergeom([a, b], c,pc^2);
    end
end
BF10 = BF10 + BF10';
logBF10 = log(BF10);
H_z_signif = (logBF10 < -3); % -3: strong evidence for H0
Hsignif = H.*~(H_z_signif);
pcorr_signif = pcorr.* ~(H_z_signif);

%% Spectral clustering (SC) on transformed graph

% Cp = reshape(diag(DCM.Cp),n,n); So, Cp is n*n (not n^2 * n^2)
Pp = 1./Cp; 
J_tf = (Pp).*(J.^2) + log(Cp) - log(pC);
J_tf (isnan(J_tf) | J_tf == Inf | J_tf ==-Inf) = 0;
J_tf_sym = (J_tf + J_tf')/2;
lap_type = 'normalized';

% SC on J_tf_sym
G = [];
param_CSC.lap_type = lap_type;  % type of Laplacian used; 'normalized' or 'combinatorial'
G.N = n;                        % number of nodes
G.W = J_tf_sym;                 % adjacanecy
G.k = 2;                        % # partitions for SC
[IDX_J(1,:)] = SC_filloutlier(G,param_CSC.lap_type,1);

% Construct block-diagonal J_pt
J_pt = zeros(n);
for i = 1:2 % for each partition
    pt = zeros(n,n);
    range = (IDX_J == i);
    pt(range,range) = J(range,range);
    J_pt = pt + J_pt;
end
%% Identigy MB nodes/states based on J partitioning

% Significantly nonzero J entries
J_nz_signif = (sign(J-1.96*sqrt(Cp)).*sign(J+1.96*sqrt(Cp))) > 0 ;
Jsignif = J.* J_nz_signif;
J_pt_signif = J_pt.*J_nz_signif;
J_btw_signif = (Jsignif - J_pt_signif);
MB_nodes_J = union(find(sum(abs(J_btw_signif),2))',find(sum(abs(J_btw_signif),1)));

%% Identify particles and their INT and MB based on J

PT1 = find(IDX_J==1); % partition 1
PT2 = find(IDX_J==2); % partition 2
MB1 = intersect(PT1,MB_nodes_J);
MB2 = intersect(PT2,MB_nodes_J);
INT1 = setdiff(PT1,MB1);
INT2 = setdiff(PT2,MB2);

%% Plot sorted J,H and pcorr

sorted = [INT1,MB1,MB2,INT2];
Jsort = Jsignif(sorted,sorted);
Hsort = Hsignif(sorted,sorted);
pcorr_sort = pcorr_signif(sorted,sorted);
sort = [{INT1},{MB1},{MB2},{INT2},{PT1},{PT2}];
Mat = [{Jsort},{Hsort},{pcorr_sort}];
Titles = [{'Effective Connectivity (J)'},...
{'Inverse Covariance (\Sigma_x^-^1)'},...
{'Partial Correlation (\rho_X_Y_|_Z)'}];

plot_HJ(Mat,sort,Titles,0,1)

%% Plot evidence for presence of MB-inducing functional connections

mask_btw_L = zeros(n,n);
mask_btw_L (INT1,PT2) = 1;
mask_btw_L (MB1,INT2) = 1;
pcorr_potentialB = pcorr(mask_btw_L==1);
logBF10_potentialB = logBF10(mask_btw_L==1);

figure;
scatter(pcorr_potentialB,logBF10_potentialB,'*')
xlabel('Partial Correlation (\rho_X_Y_|_Z)')
ylabel('LogBF10');
hold on
plot(pcorr_potentialB,-3*ones(size(pcorr_potentialB)),'-r',...
    'linewidth',2)
title('Evidence for MB-inducing functional connections')
ax = gca;ax.FontWeight = 'bold';grid on
saveas(gca,'Evidence for MB-inducing FC.png')
