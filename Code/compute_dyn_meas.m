
%% Compute dynamical measures for partitions
% dissipative timescale, solenoidal freq & kinetic energy
% Difinitions in Friston,2021: Parcels and particles: Markov blankets in the brain.

%% Identify MB states based on J partitioning

% Significantly nonzero entries of J
J_nz_signif = (sign(J-1.96*sqrt(Cp)).*sign(J+1.96*sqrt(Cp))) > 0 ;
Jsignif = J.* J_nz_signif;
J_pt_signif = J_pt.*J_nz_signif;
J_btw_signif = (Jsignif - J_pt_signif);
MB_nodes_J = union(find(sum(abs(J_btw_signif),2))',find(sum(abs(J_btw_signif),1)));

%% Identify particles and their INT and MB states based on J

PT1 = find(IDX_J==1); % partition 1
PT2 = find(IDX_J==2); % partition 2
MB1 = intersect(PT1,MB_nodes_J);
MB2 = intersect(PT2,MB_nodes_J);
INT1 = setdiff(PT1,MB1);
INT2 = setdiff(PT2,MB2);
MBs = {MB1,MB2};
tau = [];
KE = [];
freq = [];
% Compute Tau, KE, freq
for m = 1:2
    mb = MBs{m};
    Jmb = J(mb,mb);
    [e,r] = eig(Jmb);
    r     = diag(r);
    [d,jj] = sort(real(r),'descend');
    if isempty(d )== 0 
        t  = max(-8,d(min(numel(d),8 + 1)));
        nn = sum(d >= t);
        s  = r(  jj(1:nn));
        v  = e(:,jj(1:nn));
        sr = (real(s)); % or get max
        si = (imag(s)); % or get max
        tau = [tau; -1./sr];
        freq = [freq; abs(si)/(2*pi)];
        KE  = [KE;(s.*s./(-2*sr))];
    end
end
