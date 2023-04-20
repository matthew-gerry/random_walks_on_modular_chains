% pdf_plots_modRW.m

% Plot the probability distribution function derived directly from the
% master equation for a modular random walk on a long but finite chain (as
% in the function pdf_direct)

% Matthew Gerry, April 2023

%%% PARAMETERS %%%
tau = 1.0; % 1/s, tunneling element
ga_av = 1.0; % 1/s, average decoherence rate
dga = 1.0; % 1/s, difference between decoherence rates

b = 0.2; % Bias
m_list = [1,2,8]; % Segment length (even segments)

% Simulation parameters
numsites = 121; % Number of sites
dt = 0.1; % s, time step
tmax = 80; % s, max time
time = 0:dt:tmax; % time array


%%% GET PROBABILITY DISTRIBUTIONS AND AVERAGED QUANTITIES %%%

% Pre-allocate memory for results
bigPDF = zeros(numsites,length(time),length(m_list));
big_n_av = zeros(length(time),length(m_list));
big_v_av = zeros(length(time),length(m_list));
big_D_av = zeros(length(time),length(m_list));

% Get results for analogous homogeneous walk (dga=0)
[PDF_hom, sites, n_av_hom, ~, D_av_hom] = pdf_direct(2,2,b,ga_av,0,tau,numsites,dt,tmax);

for ii=1:length(m_list)
    m = m_list(ii);

    [PDF, ~, n_av, ~, D_av] = pdf_direct(m,m,b,ga_av,dga,tau,numsites,dt,tmax);

    bigPDF(:,:,ii) = PDF;
    big_n_av(:,ii) = n_av;
    big_D_av(:,ii) = D_av;
end


%% % PLOT PROBABILITY DISTRIBUTIONS %%%

% snapshot_times = floor((tmax/dt)*[0.25, 0.6, 0.99]); % Specific indices of time array at which to show PDF
snapshot_times = [200, 450, 800];

colourlist = ["#0072BD", "#D95319", "#77AC30","#7E2F8E"];
lettlist = ["(a) ", "(b) ", "(c) ", "(d) "];

figure(1)

% Plot the PDF for the homogeneous walk
subplot(2,2,1); box on; hold on
for jj=1:length(snapshot_times)
    t_snap = snapshot_times(jj);
    bar(sites, PDF_hom(:,t_snap), facecolor=colourlist(jj), facealpha=0.6, EdgeColor="none",DisplayName=strcat("$t=\;$",num2str(t_snap),"$\;s$"));
end % jj
ylim([0, 1.2*max(PDF_hom(:,snapshot_times(1)))])
yl = ylim; xl = xlim;
text(0.9*xl(1) + 0.1*xl(2), 0.1*yl(1) + 0.9*yl(2), "(a) Homogeneous", FontSize=14, Interpreter="latex")
ylabel("$P_n$",Interpreter="latex")
legend(Interpreter="latex", Location="southwest")
set(gca, fontsize=14)
hold off

% Plot the PDF for the modular walks
for ii=1:length(m_list)
    subplot(2,2,ii+1); box on; hold on
    for jj=1:length(snapshot_times)
        t_snap = snapshot_times(jj);
        bar(sites, bigPDF(:,t_snap,ii), facecolor=colourlist(jj), facealpha=0.6, EdgeColor="none");
    end % jj
    ylim([0, 1.2*max(PDF_hom(:,snapshot_times(1)))])
    yl = ylim; xl = xlim;
    if ii>1
        xlabel("$n$",Interpreter="latex")
        if rem(ii,2)==0; ylabel("$P_n$",Interpreter="latex"); end
    end
    text(0.9*xl(1) + 0.1*xl(2), 0.1*yl(1) + 0.9*yl(2), strcat(lettlist(ii+1),"$m=\;$",num2str(m_list(ii))), FontSize=14, Interpreter="latex")
    set(gca, fontsize=14)
    hold off
end % ii


%% % PLOT CUMULANTS OVER TIME %%%

figure(2)

% Mean
subplot(1,2,1); box on; hold on;
plot(time, n_av_hom./time)
for ii=1:length(m_list)
    plot(time,big_n_av(:,ii)./time')
end % ii
hold off

% Diffusion coefficient
subplot(1,2,2); box on; hold on;
plot(time, D_av_hom)
for ii=1:length(m_list)
    plot(time,big_D_av(:,ii))
end % ii
hold off