% modular_RW_CGF.m

% THIS SCRIPT IS NOT WORKING AS DESIRED - CALLS diffusionCGF_symb FUNCTION
% WHICH COULD NOT PROPERLY GET AN EXPRESSION FOR THE CGF AND THUS HAS SINCE
% BEEN MODIFIED AND RENAMED

% Investigating the properties of the cumulant generating function for the
% infinite (equivalent cycle) modular random walk. Particularly we are
% interested in the impact of bias and block_length on the CGF and
% associated cumulants.

% Matthew Gerry, February 2023

%%% CASE 1: NO BIAS

% Numerical values for parameters
tau_val = 1;
ga_av_val = 10;
dga_val = 5;
b_val = 1; % no bias

n_max = 3;
n_list = 1:n_max; % List of blocks lengths (fix two regions equal in size)
chi_list = 0:0.01:5; % For plotting the CGF

% Allocate matrices to hold values of CGF
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
[tau, ga_av, dga, b, chi, G] = diffusionCGF_symb(n,n);

J = subs(diff(G,chi),chi,0);
S = subs(diff(G,chi,2),chi,0);

% SO FAR, NOT WORKING
