function [dF,sE,sC] = spm_log_evidence_naive(qE,rE,Cp,pC)

% Return the log-evidence of a reduced model 
% (under Laplace approximation and assumption of statistical indep of connections)

% qE      - posterior expectation of full model
% rE      - retained connections of the full model = prior of reduced model
% Cp      - posterior variance of the full model (n*n, not n^2 * n^2)
% pC      - prior variance of full and reduced model (n*n, not n^2 * n^2)

% F        - reduced log-evidence: ln p(y|reduced model) - ln p(y|full model)
% [sE,sC]  - posterior expectation and covariance of reduced model
 
% Tahereh Zarghami, 2025

% Compute reduced log-evidence
%==========================================================================
% vectorize the matrices
qE  = spm_vec(qE);
qE_pt  = spm_vec(rE);
sE  = qE_pt; % posterior mean of the reduced model
pC  = spm_vec(pC);
Cp  = spm_vec(Cp);

qE_C_raw = qE - qE_pt; % removed connections (from full model)

ind = (abs(pC)>0) & (abs(qE_C_raw)>0);
ind = spm_vec(ind);

qE_C  = qE_C_raw (ind); % posterior mean of removed connections ONLY
Cp_C  = Cp(ind); % posterior variance of removed connections
Pp_C  = 1./Cp_C; % posterior precision of removed connections
pC_C  = pC(ind); % prior variance of removed connections

N   = size(qE,1);
n   = sqrt(N);

% log-evidence
%--------------------------------------------------------------------------

dF   = -0.5*( (qE_C.^2)'* Pp_C + sum(log(Cp_C)-log(pC_C)) ); % deltaF=FR-Ffull

% Posterior of the reduced model (i.e. retained connections)
%--------------------------------------------------------------------------

sE  = reshape(sE,n,n); % posterior expectation
sC  = Cp; 
sC(ind) = 0; % posterior variance
sC = reshape(sC,n,n);

