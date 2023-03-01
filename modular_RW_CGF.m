% modular_RW_CGF.m
% Numerically plotting the CGF with respect to various quantities, with and
% without bias. This is for the modular random walk on an infinite chain.

% Matthew Gerry, March 2023

% GLOBAL PARAMETERS
tau = 1;
ga_av = 10;

% HOMOGENOUS RW (SYMMETRIC/UNBIASED) - ANALYTICS
k_av = tau^2/ga_av; % classical transition rates based on model of Cao (NJP 15, 085010, 2013)
% Mean flux at steady state is trivially zero
S_av = 2*k_av; % Classic result for symmetric random walk


%%% CASE 1: plot the CGF wihtout bias at various values of delta_gamma

n = 3; % Block length for plotting against dga

% Set parameters for counting field
dchi = 0.01;
chisteps = 141; % Must be odd

dga_list = 0:4:16; % Different dga values to plot at

CGF_plot_dga = zeros(length(dga_list), chisteps); % Pre-allocate matrix of CGF values

for ii=1:length(dga_list)
    dga = dga_list(ii);

    % Chi-dressed rate matrix for the process
    [Lchi,~,~,chi] = diffusionLchi(n, n, 1, ga_av, dga, tau, dchi, chisteps);

    % Cumulant generating function
    G = CGFclassical(Lchi);
    CGF_plot_dga(ii,:) = G; % Record values
end % ii


figure; box on; hold on
plot(chi, -0.5*S_av*chi.^2, '--k',DisplayName="$\Delta\gamma=0$ (analytic)")
for ii=1:length(dga_list)
    plot(chi, real(CGF_plot_dga(ii,:)),DisplayName=strcat("$\Delta\gamma$ = ",num2str(dga_list(ii))))
end % ii
legend(Interpreter="latex",Location="south")
xlim([min(chi),max(chi)])
ylim([-0.04,0.005])
xlabel("$\chi$",Interpreter="latex")
ylabel("$\mathcal{G}(\chi)$",Interpreter="latex")
title(strcat("$n_A=n_B=\;$",num2str(n),", $\bar{\gamma} =\;$",num2str(ga_av)),Interpreter="latex")
set(gca, fontsize=14)
hold off


%%% CASE 2: plot the CGF for various values of the bias ratio b, with fixed delta gamma

b_list = 1:1:4; % Various values of the ratio of forward to reverse rates
dga_b = 5; % dga value for plotting at varying b

CGF_plot_b = zeros(length(b_list), chisteps); % Pre-allocate matrix of CGF values

for ii=1:length(b_list)
    b = b_list(ii);

    % Chi-dressed rate matrix for the process
    [Lchi,~,~,chi] = diffusionLchi(n, n, b, ga_av, dga_b, tau, dchi, chisteps);

    % Cumulant generating function
    G = CGFclassical(Lchi);
    CGF_plot_b(ii,:) = G; % Record values
end % ii

figure;
subplot(1,2,1); box on; hold on
for ii=1:length(b_list)
    plot(chi, imag(CGF_plot_b(ii,:)),DisplayName=strcat("$b$ = ",num2str(b_list(ii))));
end % ii
legend(Interpreter="latex",Location="north",Orientation="horizontal")
xlim([min(chi),max(chi)])
xlabel("$\chi$",Interpreter="latex")
ylabel("Im$[\mathcal{G}(\chi)]$",Interpreter="latex")
% title(strcat("$n_A=n_B=\;$",num2str(n),", $\bar{\gamma} =\;$",num2str(ga_av)),Interpreter="latex")
set(gca, fontsize=14)
hold off

subplot(1,2,2); box on; hold on
for ii=1:length(b_list)
    plot(chi, real(CGF_plot_b(ii,:)),DisplayName=strcat("$b$ = ",num2str(b_list(ii))));
end % ii
% legend(Interpreter="latex")
xlim([min(chi),max(chi)])
xlabel("$\chi$",Interpreter="latex")
ylabel("Re$[\mathcal{G}(\chi)]$",Interpreter="latex")
title(strcat("$n_A=n_B=\;$",num2str(n),", $\bar{\gamma} =\;$",num2str(ga_av),", $\Delta\gamma=\;$",num2str(dga_b)),Interpreter="latex")
set(gca, fontsize=14)
hold off
