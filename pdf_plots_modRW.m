% pdf_plots_modRW.m

% Plot the probability distribution function derived directly from the
% master equation for a modular random walk on a long but finite chain (as
% in the function pdf_direct)

% Matthew Gerry, April 2023

%%% PARAMETERS %%%
tau = 1; % 1/s, tunneling element
ga_av = 1.0; % 1/s, average decoherence rate
dga = 1.0; % 1/s, difference between decoherence rates

b = 0.2; % Bias
m_list = [1,2,4,8]; % Segment length (even segments)

% Simulation parameters
numsites = 161; % Number of sites
dt = 0.2; % s, time step
tmax = 120; % s, max time
time = 0:dt:tmax; % time array


%%% GET PROBABILITY DISTRIBUTIONS AND AVERAGED QUANTITIES %%%

% Pre-allocate memory for results
bigPDF = zeros(numsites,length(time),length(m_list));
big_n_av = zeros(length(time),length(m_list));
big_v_av = zeros(length(time),length(m_list));
big_D_av = zeros(length(time),length(m_list));
big_C3 = zeros(length(time),length(m_list));
big_C4 = zeros(length(time),length(m_list));

% Get results for analogous homogeneous walk (dga=0)
[PDF_hom, sites, n_av_hom, ~, D_av_hom, C3_hom, C4_hom] = pdf_direct(2,2,b,ga_av,0.0,tau,numsites,dt,tmax);

for ii=1:length(m_list)
    m = m_list(ii);

    % Derive probability distributions and cumulants using PDF_direct function
    [PDF, ~, n_av, ~, D_av, C3, C4] = pdf_direct(m,m,b,ga_av,dga,tau,numsites,dt,tmax);

    bigPDF(:,:,ii) = PDF;
    big_n_av(:,ii) = n_av;
    big_D_av(:,ii) = D_av;
    big_C3(:,ii) = C3;
    big_C4(:,ii) = C4;
end

% Define some colours, linestyles, and labels for plotting
colourlist = ["#0072BD", "#D95319", "#77AC30","#7E2F8E"];
lettlist = ["(a) ", "(b) ", "(c) ", "(d) "];
ls_list = ["-","--",":","-."];


%%% PLOT PROBABILITY DISTRIBUTIONS %%%

snapshot_times = [50, 175, 400];  % Specific indices of time array at which to show PDF (arb.)

figure(1)

% Plot the PDF for the modular walks
for ii=1:length(m_list)
    subplot(2,2,ii); box on; hold on
    for jj=1:length(snapshot_times)
        t_snap = snapshot_times(jj);
        hom_curve = plot(sites, PDF_hom(:,t_snap), '-k', linewidth=0.2);
        bar(sites, bigPDF(:,t_snap,ii), facecolor=colourlist(jj), facealpha=0.6, EdgeColor="none",DisplayName=strcat("$t=\;$",num2str(dt*t_snap)));
        hom_curve.Annotation.LegendInformation.IconDisplayStyle = "off";
    end % jj
%     xlim([-10,sites(end)]) % If high bias
    ylim([0, 1.2*max(bigPDF(:,snapshot_times(1),ii))])
    yl = ylim; xl = xlim;
    if ii>2; xlabel("$n$",Interpreter="latex"); end
    if rem(ii,2)~=0; ylabel("$P_n$",Interpreter="latex"); end
    if ii==1; legend(Interpreter="latex", Location="southwest"); end
    text(0.9*xl(1) + 0.1*xl(2), 0.1*yl(1) + 0.9*yl(2), strcat(lettlist(ii),"$m=\;$",num2str(m_list(ii))), FontSize=14, Interpreter="latex")
    set(gca, fontsize=14)
    hold off
end % ii


%%% PLOT CUMULANTS OVER TIME %%%

figure(2)

% Mean
subplot(2,2,1); box on; hold on;
plot(time, n_av_hom./time, '--k', DisplayName="Homogeneous")
for ii=1:length(m_list)
    plot(time,big_n_av(:,ii)./time', ls_list(ii), color=colourlist(ii), linewidth=1.5, DisplayName=strcat("$m=\;$",num2str(m_list(ii))))
end % ii
xlim([0, 100])
ylabel("$\mathcal{C}_1$", Interpreter="latex")
legend(Interpreter="latex", Location="southeast")
set(gca,fontsize=14)
hold off

% Diffusion coefficient times 2 (scaled variance)
subplot(2,2,2); box on; hold on;
plot(time, 2*D_av_hom, '--k')
for ii=1:length(m_list)
    plot(time,2*big_D_av(:,ii), ls_list(ii), color=colourlist(ii), linewidth=1.5)
end % ii
xlim([0, 100])
ylabel("$\mathcal{C}_2$", Interpreter="latex")
set(gca,fontsize=14)
hold off

% Scaled skewness
subplot(2,2,3); box on; hold on;
plot(time, C3_hom, '--k')
for ii=1:length(m_list)
    plot(time,big_C3(:,ii), ls_list(ii), color=colourlist(ii), linewidth=1.5)
end % ii
xlim([0, 100])
ylabel("$\mathcal{C}_3$", Interpreter="latex")
xlabel("$t$", Interpreter="latex")
set(gca,fontsize=14)
hold off

% Scaled kurtosis
subplot(2,2,4); box on; hold on;
plot(time, C4_hom, '--k')
for ii=1:length(m_list)
    plot(time,big_C4(:,ii), ls_list(ii),color=colourlist(ii), linewidth=1.5)
end % ii
xlim([0, 100])
ylabel("$\mathcal{C}_4$", Interpreter="latex")
xlabel("$t$", Interpreter="latex")
set(gca,fontsize=14)
hold off