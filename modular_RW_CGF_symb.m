% modular_RW_CGF.m

% Investigating the properties of the cumulant generating function for the
% infinite (equivalent cycle) modular random walk. Particularly we are
% interested in the impact of bias and block_length on the CGF and
% associated cumulants.

% We expand the CGF in chi, and then expand each of the coefficient in
% powers of the variables of interest: delta_gamma and b.

% Matthew Gerry, March 2023

%%% CASE 1: NO BIAS

% Numerical values for parameters
% tau_val = 1;
% ga_av_val = 10;
% dga_val = 5;
% b_val = 1; % no bias
% 
% n_max = 3;
% n_list = 1:n_max; % List of blocks lengths (fix two regions equal in size)
% chi_list = 0:0.01:5; % For plotting the CGF
% 
% % Allocate matrices to hold values of CGF
% G_re_vals = zeros(n_max, length(chi_list));
% G_im_vals = zeros(n_max, length(chi_list));

%%% JUST PLOT THE CGF NUMERICALLY, THERE IS NO NEED TO USE SYMBOLIC MANIPULATION

% syms tau ga_av dga b chi

% for n=n_list
%     [tau, ga_av, dga, b, chi, G] = diffusionCGF_symb(n,n);
% 
%     % Substitute numerical values
%     G_plot = subs(G, [tau, ga_av, dga, b], [tau_val, ga_av_val, dga_val, b_val]);
%     G_plot = double(subs(G_plot,chi,chi_list));
%     
%     % Record values
%     G_re_vals(n,:) = real(G_plot);
%     G_im_vals(n,:) = imag(G_plot);
% end % n
%     
% figure
% subplot(1,2,1); hold on; box on
% for n=n_list
%     plot(chi_list, G_im_vals(n,:))
% end
% hold off
% 
% subplot(1,2,2); hold on; box on
% for n=n_list
%     plot(chi_list, G_re_vals(n,:))
% end
% hold off


%%% ANALYTIC EXPRESSIONS FOR CUMULANTS, WITH BIAS
n = 2; % For now
[tau, ga_av, dga, b, chi, Lchi] = diffusionLchi_symb(n,n);

Lchi_eigs = eig(Lchi); % Eigenvalues of the chi-dressed rate matrix

% Verify which eigenvalue is the scaled CGF by plugging in chi=0 and
% identifying the one that vanishes - note there is discontinuous behavour
% here at zero bias - it appears to be the first element of Lchi_eigs for
% arbitrarily small but finite bias, but then its the second element at
% exactly zero bias
Lchi_mixed = simplify(subs(Lchi, [tau, ga_av, dga, b, chi], [1, 10, 5, 0.00001, 0]));
mixed_eigs = double(eig(Lchi_mixed));
CGF_index = find(mixed_eigs==0);

CGF = Lchi_eigs(CGF_index);

% Let's begin the expansions
% Mean current
J = subs(diff(CGF, chi), chi, 0); % Symbolic function for the current




% SO FAR, NOT WORKING


