function [H,R] = spm_ness_mod(J,G)
% This is a MODIFIED version of spm_ness.m in SPM12
% Modified by Tahereh Zarghami, 2025

% Evaluation of hessian and solenoidal operators at NESS
% FORMAT [H,R]     = spm_ness(J,G)
% FORMAT [H,R,J,G] = spm_ness(J,G)
% J  - Jacobian (dfdx)
% G  - diffusion tensor (amplitude of random fluctuations)
%
% H  - Hessian matrix (i.e., precision of a Gaussian density)
% R  - Skew symmetric solenoidal operator (-Q')
%
% if called with four output arguments, complex forms are returned
%__________________________________________________________________________
% This routine evaluates the Hessian (i.e., precision) of a nonequilibrium
% steady-state density (using a local linear approximation, under Gaussian
% assumptions). This is evaluated  under linear constraints on the
% solenoidal term of a Helmholtz decomposition. In short, given the flow
% (encoded by the systems Jacobian) and amplitude of random fluctuations,
% one can evaluate the steady-state density under nonequilibrium dynamics
% implied by solenoidal flow.
%
% There are additional notes using symbolic maths and numerical examples in
% the main body of the script.
%
% flow constraints (Jacobian J)(R = -Q')
%--------------------------------------------------------------------------
% where flow   f = (R + G)*d log(p(x))/dx and
% log(p(x))      = -(1/2)*x'*H*x =>
% d log(p(x))/dx = -H*x =>
% df/dx = J      = -(R + G)*H =>
% H              = -(R + G)\J =>
% J*R + R*J'     = J*G - G*J'
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ness.m 7910 2020-07-28 19:16:17Z karl $


% solve for solenoidal (R) operator J*R + R*J' = J*G - G*J'
%==========================================================================

% Original in spm:
%     n  = size(J,1);
%     I  = eye(n,n);
%     X  = kron(I,J) + kron(conj(J),I);
%     Y  = spm_vec(J*G - G*J');
%     R  = reshape(X\Y,n,n);

% Modified: solve sylvester eq using Matlab's built-in function:
Y  = (J*G - G*J');
R = sylvester(J,J',Y);

% precision (inverse covariance) of steady-state density
%----------------------------------------------------------------------
H  = -(R + G)\J;
    
