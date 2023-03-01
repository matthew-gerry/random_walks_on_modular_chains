% modular_RW_CGF.m
% Numerically plotting the CGF with respect to various quantities, with and
% without bias. This is for the modular random walk on an infinite chain.

% Matthew Gerry, March 2023

%%% Plot the CGF wihtout bias at various values of delta_gamma

n_dga = 3; % Block length for plotting against dga

% Parameters
tau = 1;
ga_av = 10;
dchi = 0.01;
chisteps = 61; % Must be odd

dga_list = 0:2:10; % Different dga values to plot at

CGF_plot = zeros(length(dga_list), chisteps); % Pre-allocate matrix of CGF values

for ii=1:length(dga_list)
    dga = dga_list(ii);

    % Chi-dressed rate matrix for the process
    [Lchi,~,~,chi] = diffusionLchi(n_dga, n_dga, 1, ga_av, dga, tau, dchi, chisteps);

    % Cumulant generating function
    G = CGFclassical(Lchi);
    CGF_plot(ii,:) = G; % Record values
end % ii

figure; box on; hold on
for ii=1:length(dga_list)
    plot(chi, real(CGF_plot(ii,:)),DisplayName=strcat("$\Delta\gamma = ",num2str(dga_list(ii))));
end % ii
legend(Interpreter="latex")
xlabel(
hold off