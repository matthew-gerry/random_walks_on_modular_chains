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
    [Lchi,~,~,chi] = diffusionLchi(n, n, 0, ga_av, dga, tau, dchi, chisteps);

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


%%% CASE 2: plot the CGF for various values of the bias lr, with fixed delta gamma

bias_list = [0, 0.5, 1, 2]; % Various values of the log-ratio of forward to reverse rates
dga_b = 5; % dga value for plotting at varying bias

CGF_plot_b = zeros(length(bias_list), chisteps); % Pre-allocate matrix of CGF values

for ii=1:length(bias_list)
    bias = bias_list(ii);

    % Chi-dressed rate matrix for the process
    [Lchi,~,~,chi] = diffusionLchi(n, n, bias, ga_av, dga_b, tau, dchi, chisteps);

    % Cumulant generating function
    G = CGFclassical(Lchi);
    CGF_plot_b(ii,:) = G; % Record values
end % ii

figure;
subplot(1,2,1); box on; hold on
for ii=1:length(bias_list)
    plot(chi, imag(CGF_plot_b(ii,:)),DisplayName=strcat("Bias: ",num2str(bias_list(ii))));
end % ii
legend(Interpreter="latex",Location="north",Orientation="horizontal")
xlim([min(chi),max(chi)])
xlabel("$\chi$",Interpreter="latex")
ylabel("Im$[\mathcal{G}(\chi)]$",Interpreter="latex")
% title(strcat("$n_A=n_B=\;$",num2str(n),", $\bar{\gamma} =\;$",num2str(ga_av)),Interpreter="latex")
set(gca, fontsize=14)
hold off

subplot(1,2,2); box on; hold on
for ii=1:length(bias_list)
    plot(chi, real(CGF_plot_b(ii,:)),DisplayName=strcat("Bias: ",num2str(bias_list(ii))));
end % ii
% legend(Interpreter="latex")
xlim([min(chi),max(chi)])
xlabel("$\chi$",Interpreter="latex")
ylabel("Re$[\mathcal{G}(\chi)]$",Interpreter="latex")
title(strcat("$n_A=n_B=\;$",num2str(n),", $\bar{\gamma} =\;$",num2str(ga_av),", $\Delta\gamma=\;$",num2str(dga_b)),Interpreter="latex")
set(gca, fontsize=14)
hold off


%% % CASE 3: compare to the analytic expression for the CGF when n=1 for both regions

bias_list_compare = [0, 0.5, 4];

% Set parameters for counting field
dchi = 0.04;
chisteps = 61; % Must be odd


CGF_num = zeros(length(bias_list_compare), chisteps); % Pre-allocate matrices of CGF values
CGF_ana = zeros(length(bias_list_compare), chisteps);
for ii=1:length(bias_list_compare)
    b = bias_list_compare(ii);

    [Lchi,~,~,chi] = diffusionLchi(1, 1, b, ga_av, dga_b, tau, dchi, chisteps);

    % Cumulant generating function determined numerically
    G = CGFclassical(Lchi);
    CGF_num(ii,:) = G; % Record values

    % Analytic expression for the CGF derived in modular_RW_CGF_symb.m
    sqrtarg = (cosh(0.5*b+1i*chi)).^2 - 0.25i*(dga_b/ga_av)^2*sin(chi).*sinh(b+1i*chi);
    G_ana = 2*tau^2*exp(-b/2)*(sqrt(sqrtarg) - cosh(b/2))/(ga_av*(1-0.25*(dga_b/ga_av)^2));

    CGF_ana(ii,:) = G_ana;

end

colour_list = ['b','r','m'];

% Numerical values are plotted with marker, analytic with dashed curves

figure
subplot(1,2,1); hold on; box on
for ii=1:length(bias_list_compare)
    numcurve = plot(chi, imag(CGF_num(ii,:)), 'x', Color=colour_list(ii));
    numcurve.Annotation.LegendInformation.IconDisplayStyle = "off";
    plot(chi, imag(CGF_ana(ii,:)), '--', Color=colour_list(ii), DisplayName=strcat("$b=\;$",num2str(bias_list_compare(ii))))
end
xlabel("$\chi$", Interpreter="latex")
ylabel("Im[$\mathcal{G}(\chi)$]",Interpreter="latex")
legend(Interpreter="latex",Location="northwest")
set(gca, fontsize=14)
hold off

subplot(1,2,2); hold on; box on
for ii=1:length(bias_list_compare)
    numcurve = plot(chi, real(CGF_num(ii,:)), 'x', Color=colour_list(ii));
    numcurve.Annotation.LegendInformation.IconDisplayStyle = "off";
    plot(chi, real(CGF_ana(ii,:)), '--', Color=colour_list(ii), DisplayName=strcat("$b=\;$",num2str(bias_list_compare(ii))))
end
xlabel("$\chi$", Interpreter="latex")
ylabel("Re[$\mathcal{G}(\chi)$]",Interpreter="latex")
% legend(Interpreter="latex")
set(gca, fontsize=14)
hold off

