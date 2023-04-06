% modular_RW_CGF.m

% Investigating the properties of the cumulant generating function for the
% infinite (equivalent cycle) modular random walk. Particularly we are
% interested in the impact of bias and block_length on the CGF and
% associated cumulants.

% We expand the CGF in chi, and then expand each of the coefficient in
% powers of the variables of interest: delta_gamma and b.

% Matthew Gerry, March 2023

%%% ANALYTIC EXPRESSIONS FOR CUMULANTS, WITH BIAS

n = 1; % For now
[tau, ga_av, dga, b, chi, Lchi] = diffusionLchi_symb(n,n);
% Lchi = subs(Lchi, [tau, ga_av], [1,10]); % Substitute some numerical values to speed up code

Lchi_eigs = eig(Lchi); % Eigenvalues of the chi-dressed rate matrix

% Verify which eigenvalue is the scaled CGF by plugging in chi=0 and
% identifying the one that vanishes - note there is discontinuous behavour
% here at zero bias - it appears to be the first element of Lchi_eigs for
% arbitrarily small but finite bias, but then its the second element at
% exactly zero bias
% Lchi_mixed = simplify(subs(Lchi, chi, 0));
evaluated_eigs = double(subs(Lchi_eigs,[tau, ga_av, dga, b, chi],[1,10,5,0.5,0]));
CGF_index = find(evaluated_eigs==0);

CGF = Lchi_eigs(CGF_index);

% Let's begin the expansions

% Mean current
J = -1i*subs(diff(CGF, chi), chi, 0); % Symbolic function for the current
% Mean current dependence on delta_gamma
J_dga_1 = diff(J,dga);

% Variance
S = -simplify(subs(diff(CGF,chi,2),chi,0)); % Symbolic function for the variance

% Skewness
C3 = simplify(1i*subs(diff(CGF,chi,3),chi,0));

% Kurtosis
C4 = simplify(subs(diff(CGF,chi,4),chi,0));


% % Plot the variance as a check
% S_plot = subs(S,[tau, ga_av, dga],[1.0,10.0,5.0]);
% figure
% fplot(S_plot,[0,10]);

%% Analytic results with bias, longer segments
% Just the steady state probability distribution

n = 2; % Won't run with more than 2
[tau, ga_av, dga, b, chi, Lchi] = diffusionLchi_symb(n,n);

L = simplify(subs(Lchi,chi,0));

[V, D] = eig(L);

D = diag(D);

ind = find(D==0);
P = V(:, ind);
P = simplify(P/sum(P));

