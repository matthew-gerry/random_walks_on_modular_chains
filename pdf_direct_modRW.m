% pdf_direct_modRW.m

% Calculate the probability distribution function directly for a modular
% random walk on a long chain, given an initial state (at site zero, by
% choice). Do this by solving the master equation numerically.

% Matthew Gerry, April 2023

%%% SET UP PROBLEM %%%

% Constant parameters
tau = 1.0; % 1/s, "tunnel coupling"
ga_av = 1.0; % 1/s, average "decoherence rate"
dga = 1.0; % 1/s, difference between decoherence rates

m = 4; % Segment length
b = 0.5; % Bias
% b  = -log((ga_av-0.5*dga)/(ga_av+0.5*dga));


numsites = 201; % Total number of sites
p0 = zeros([numsites,1]); p0((numsites+1)/2) = 1; % Initial state

% Calculate the rate matrix for this random walk
[L, n, block_types] = L_explicit(m,m,b,ga_av,dga,tau,numsites);

%%% SOLVE MASTER EQUATION NUMERICALLY %%%
dt = 0.1;
tmax = 100;
time = 0:dt:tmax;

p = zeros(numsites,length(time)); % Pre-allocate time-series of prob dist

for ii=1:length(time)
    t = time(ii);
    p(:,ii) = expm(L*t)*p0;
end

% Prob dist plots
figure; hold on
plot(n, p(:,20))
plot(n, p(:,200))
plot(n, p(:,400))
plot(n, p(:,800))
hold off

%%% STATISTICS OF n %%% 
dpdt = L*p;

n_av = sum(p.*repmat(n',[1,length(time)]));
v_av = sum(dpdt.*repmat(n',[1,length(time)]));
D_av = sum(dpdt.*repmat((n.^2)',[1,length(time)])) - 2*n_av.*v_av;

figure
plot(time, n_av)

figure
plot(time, v_av)

figure
plot(time, D_av)
