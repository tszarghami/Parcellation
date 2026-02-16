
function [IDX_SC, Uk, Dk, time_SC] = SC_filloutlier(G,lap_type,outlier_flag)
% This is a MODIFIED version of SC.m code of Nicolas Tremblay & Gilles Puy (CSCbox)
% SC.m is part of the CSCbox (Compressive Spectral Clustering toolbox)
% Copyright (C) 2016 Nicolas Tremblay, Gilles Puy.

if strcmp(lap_type,'normalized')
    normBOOL=true;
    G.lap_type='normalized';
    G.L = create_laplacian(G.W,G.lap_type);
elseif strcmp(lap_type,'combinatorial')
    normBOOL=false;
    G.lap_type='combinatorial';
    G.L = create_laplacian(G.W,G.lap_type);
else
    error(['SC.m : lap_type should be ''normalized'' or ''combinatorial''']);
end

% eigs
if isfield(G, 'G.Dk') % if it was already calculated, do not do it again
    Dk=G.Dk;
    Uk=G.Uk;
    time_SC.eigs=G.time_SC_eigs;
else % run eigs
    tic;
    opts.isreal=1;%opts.issym=1;opts.maxit=10000;
    [Uk,Dk]=eigs(G.L,G.k,'SA',opts);Dk=diag(Dk);
    time_SC.eigs=toc;
    % This is the modified part:
    if outlier_flag
    % Refs: Depavia2020:Info hidden in Fiedler vector & Tam2020:Fiedler Regularization:
    [Uk(:,2),TF] = filloutliers(Uk(:,2),'clip','median');
    N_Fiedler_outliers = sum(TF);
    end
end

if normBOOL % kmeans on Yk
    tic;
    Yk=Uk./repmat(sqrt(sum(Uk.^2,2)),1,G.k);
    IDX_SC = kmeans(Yk, G.k,'Replicates',20);
    time_SC.kmeans=toc;
else % kmeans on Uk
    tic;
    IDX_SC = kmeans(Uk, G.k,'Replicates',20);
    time_SC.kmeans=toc;
end

time_SC.total=time_SC.eigs+time_SC.kmeans;
