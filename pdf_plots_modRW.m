% pdf_plots_modRW.m

% Plot the probability distribution function derived directly from the
% master equation for a modular random walk on a long but finite chain (as
% in the function pdf_direct)

% Matthew Gerry, April 2023

%%% PARAMETERS %%%
tau = 1.0; % 1/s, tunneling element
ga_av = 1.0; % 1/s, average decoherence rate
dga = 1.0; % 1/s, difference between decoherence rates

b = 0.5; % Bias
m_list = [1,2,4]; % Segment length (even segments)

% Simulation parameters
numsites = 151; % Number of sites
dt = 0.1; % s, time step
tmax = 60; % s, max time
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


%%% PLOT PROBABILITY DISTRIBUTIONS %%%

snapshot_times = floor((tmax/dt)*[0.1, 0.5, 0.99]); % Specific indices of time array at which to show PDF

figure(1)

% Plot the PDF for the homogeneous walk
subplot(1,length(m_list)+1,1); box on; hold on
for jj=snapshot_times
    plot(sites, PDF_hom(:,jj));
end % jj

% Plot the PDF for the modular walks
for ii=1:length(m_list)
    subplot(1,length(m_list)+1,ii+1); box on; hold on
    for jj=snapshot_times
        plot(sites, bigPDF(:,jj,ii));
    end % jj
    hold off
end % ii


%%% PLOT CUMULANTS OVER TIME %%%

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