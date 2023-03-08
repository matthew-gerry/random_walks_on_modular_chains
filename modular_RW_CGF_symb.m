% modular_RW_CGF.m

% Investigating the properties of the cumulant generating function for the
% infinite (equivalent cycle) modular random walk. Particularly we are
% interested in the impact of bias and block_length on the CGF and
% associated cumulants.

% We expand the CGF in chi, and then expand each of the coefficient in
% powers of the variables of interest: delta_gamma and b.

% Matthew Gerry, March 2023

%%% ANALYTIC EXPRESSIONS FOR CUMULANTS, WITH BIAS

n = 2; % For now
[tau, ga_av, dga, b, chi, Lchi] = diffusionLchi_symb(n,n);
Lchi = subs(Lchi, [tau, ga_av], [1,10]); % Substitute some numerical values to speed up code

Lchi_eigs = eig(Lchi); % Eigenvalues of the chi-dressed rate matrix

% Verify which eigenvalue is the scaled CGF by plugging in chi=0 and
% identifying the one that vanishes - note there is discontinuous behavour
% here at zero bias - it appears to be the first element of Lchi_eigs for
% arbitrarily small but finite bias, but then its the second element at
% exactly zero bias
Lchi_mixed = simplify(subs(Lchi, chi, 0));
mixed_eigs = eig(Lchi_mixed);
CGF_index = find(mixed_eigs==0);

CGF = Lchi_eigs(CGF_index);

% Let's begin the expansions
% Mean current
J = -1i*subs(diff(CGF, chi), chi, 0); % Symbolic function for the current

% Mean current dependence on delta_gamma
J_dga_1 = diff(J,dga);
J_dga_2 = diff(J_dga_1, dga);

% SO FAR, TOO SLOW


