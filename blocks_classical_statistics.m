%% Investigate the statistics of a random walk with rates varying between
% "blocks" of sites (classical analogue of the system studied by Francisco)

% Start with a symmetric random walk, alternating two-site-long blocks of
% slower and faster rates r and s=a*r

% We technically examine an equivalent system which is a closed loop of 
% four states, in which we count the number of steps taken in the clockwise
% direction.

% Symbolic manipulation of the chi-dressed master equation
syms W
syms r a chi

% W = r*[[-2, exp(-1i*chi), 0, exp(1i*chi)];
%        [exp(1i*chi), -1-a, a*exp(-1i*chi), 0];
%        [0, a*exp(1i*chi), -2*a, a*exp(-1i*chi)];
%        [exp(-1i*chi), 0, a*exp(1i*chi), -1-a]];
W = r*[[-2, exp(-4i*chi), 0, 1];
       [exp(4i*chi), -1-a, a, 0];
       [0, a, -2*a, a];
       [1, 0, a, -1-a]];

[V,D] = eig(W); % Eigenvalues and eigenvectors

% Steady state - not sure if this makes sense given the model
syms G % Scaled cumulant generating function at steady state
d = diag(D);
% This method of identifying the dominant eigenvalue is probably slowing
% down the code a lot
CGFindex = find(real(double(subs(d,[r,a,chi],[1,0.5,0])))==0); % Arb. values subbed for r and a
G = d(CGFindex);

J = -1i*subs(diff(G,chi),chi,0); % Mean current at steady state
S = -subs(diff(G,chi,2),chi,0); % Variance in current at steady state

% Plug in some values and plot
r_val = 1; % arb
a_list = 0.01:0.3:20;
% The program throws an error if you try to plot at a=1 (I think there is a
% divergence due to extra symmetries in the Liouvillian) - so just offset
% the a_list so as not to include a=1. Similar with a=0.

S_plot = subs(S,r,r_val);
S_plot = double(subs(S_plot,a,a_list));

S_analytic = 4*a_list*r_val^2./(r_val*(1+a_list));

% Now do the analysis for the homogeneous random walk whose transition rate
% between all states is equal to the average of the two in the
% heterogeneous model. This should also be possible to calculate 
% analytically without counting (just a simple random walk).

% % This problem simplifies to the case of a two-state network
% 
% syms Wav chif chib
% 
% Wav = 0.5*r*(1+a)*[[-1, exp(1i*chif) + exp(1i*chib)];
%                    [exp(1i*chif) + exp(1i*chib), -1]];
% 
% [Vav, Dav] = eig(Wav);
% 
% syms Gav
% Gav = -Dav(1,1);
% 
% Javf = -1i*subs(diff(Gav,chif),[chif,chib],[0,0]);
% Javb = -1i*subs(diff(Gav,chib),[chif,chib],[0,0]);
% Jav_product = -subs(diff(diff(Gav,chif),chib),[chif,chib],[0,0]);
% 
% Savf = -subs(diff(Gav,chif,2),[chif,chib],[0,0]);
% Savb = -subs(diff(Gav,chib,2),[chif,chib],[0,0]);
% 
% Sav = Savf + Savb - 2*Jav_product;
% Sav_plot = subs(Sav,r,r_val);
% Sav_plot = double(subs(Sav_plot,a,a_list));


%
% Of course, for a continuous-time, symmetric, homogeneous random walk, the
% variance in the current is simply given by the transition rate (more
% precisely, it is 2 times the diffusion coefficient, divide out the
% spatial scale since we are not worried about continuous space here, which
% is simply 1/tau, tau the characteristic time between transitions)
Sav_analytic = 2*0.5*r_val*(1+a_list);

% Also plot the case where the rates are both equal to the sum (as a check)
% Ssum_analytic = r_val*(1+a_list);


% Plotting
figure; hold on; box on
plot(a_list, S_plot, '.b',"DisplayName", "Alternating blocks, $r_1$ and $r_2$")
plot(a_list, S_analytic, '-m', "DisplayName", "Alternating blocks, analytic expression")
% plot(a_list, Sav_plot)
plot(a_list, Sav_analytic, '--r', 'DisplayName', "Homogeneous at average rate $(r_1+r_2)/2$")
% plot(a_list, Ssum_analytic, 'DisplayName', "Homogeneous at sum of rates $(r_1+r_2)$")
xlabel("$r_2/r_1$","Interpreter","latex")
ylabel("$\langle\langle J^2\rangle\rangle$","Interpreter","latex")
legend('Interpreter','latex', Location='northwest')
set(gca,'Fontsize', 14)

hold off

%% Now the asymmetric case
% r the forward rate for one block, a*r the forwards rate for the other
% b the ratio of backwards to forwards rates, both blocks

% Symbolic manipulation of the chi-dressed master equation
syms W2
syms r a b chi

% W2 = r*[[-1-b, b*exp(-1i*chi), 0, exp(1i*chi)];
%        [exp(1i*chi), -a-b, a*b*exp(-1i*chi), 0];
%        [0, a*exp(1i*chi), -a*(1+b), a*b*exp(-1i*chi)];
%        [b*exp(-1i*chi), 0, a*exp(1i*chi), -1-a*b]];

W2 = r*[[-1-b, b*exp(-1i*chi), 0, 1];
       [exp(1i*chi), -a-b, a*b, 0];
       [0, a, -a-a*b, a*b];
       [b, 0, a, -1-a*b]];

[V2,D2] = eig(W2); % Eigenvalues and eigenvectors

% Steady state - not sure if this makes sense given the model
syms G2 % Scaled cumulant generating function at steady state
G2 = D2(1,1); % Just hard-coding that it's the first element (based on checks)
          % Difficult to soft code - this could break the code though

J2 = -1i*subs(diff(G2,chi),chi,0); % Mean current at steady state
S2 = -subs(diff(G2,chi,2),chi,0); % Variance in current at steady state

% Plug in some values and plot
r_val = 1; % arb
b_val = 0.1; % Forward rates twice as fast as reverse rates
a_list = 0.01:0.15:20;

J2_plot = subs(J2,[r,b],[r_val,b_val]);
J2_plot = double(subs(J2_plot,a,a_list));
S2_plot = subs(S2,[r,b],[r_val,b_val]);
S2_plot = double(subs(S2_plot,a,a_list));

% Now do the analysis for the homogeneous random walk whose forward (reverse) 
% transition rate between all states is equal to the average of the two in
% the heterogeneous model. This should also be possible to calculate 
% analytically without counting (just a simple random walk).
% This problem simplifies to the case of a three-state network

syms Wav2

% Wav2 = 0.5*r*(1+a)*[[-1-b, exp(1i*chi) + b*exp(-1i*chi)];
%                     [exp(1i*chi) + b*exp(-1i*chi), -1-b]];

Wav2 = 0.5*r*(1+a)*[[-1-b, b*exp(-1i*chi), 1];
                    [exp(1i*chi), -1-b, b];
                    [b, 1, -1-b]];

[Vav2, Dav2] = eig(Wav2);

syms Gav2
Gav2 = Dav2(1,1);

Jav2 = -real(1i*subs(diff(Gav2,chi),chi,0));
Sav2 = -subs(diff(Gav2,chi,2),chi,0);


Jav2_plot = subs(Jav2,[r,b],[r_val,b_val]);
Jav2_plot = double(subs(Jav2_plot,a,a_list));
Sav2_plot = subs(Sav2,[r,b],[r_val,b_val]);
Sav2_plot = double(subs(Sav2_plot,a,a_list));


% Plotting
figure
subplot(121); hold on; box on
plot(a_list, J2_plot, "DisplayName", "Alternating blocks, $r_1$ and $r_2$")
plot(a_list, Jav2_plot, 'DisplayName', "Homogeneous at average rate $(r_1+r_2)/2$")
xlabel("$r_2/r_1$","Interpreter","latex")
ylabel("$\langle J\rangle$","Interpreter","latex")
legend('Interpreter','latex', Location='northwest', Orientation='horizontal')
set(gca,'Fontsize', 14)
hold off

subplot(122); hold on; box on
plot(a_list, S2_plot, "DisplayName", "Alternating blocks, $r_1$ and $r_2$")
plot(a_list, Sav2_plot, 'DisplayName', "Homogeneous at average rate $(r_1+r_2)/2$")
xlabel("$r_2/r_1$","Interpreter","latex")
ylabel("$\langle\langle J^2\rangle\rangle$","Interpreter","latex")
% legend('Interpreter','latex', Location='northwest')
set(gca,'Fontsize', 14)
hold off

%% Check: calculate mean flux (nonequilibrium) without counting (from steady state)
% Try instead keeping the average rate between the two regions fixed

% Symbolic manipulation with the rate matrix
syms L
syms r d b

% r: the average forwards rate between the two regions
% d: the difference in forwards rates between the two regions
% b: the ratio of backwards rate to forwards rate (common between the two regions)
 
% L = r*[[-1-b, b, 0, 1];
%        [1, -a-b, a*b, 0];
%        [0, a, -a*(1+b), a*b];
%        [b, 0, a, -1-a*b]];

L = [[-(1+b)*r + (1+b)*d/2, b*(r-d/2),            0,                       r-d/2];
     [r-d/2,                -(1+b)*r - (1-b)*d/2, b*(r+d/2),               0];
     [0,                    r+d/2,                -(1+b)*r - (1+b)*d/2 ,   b*(r+d/2)];
     [b*(r-d/2),            0,                    r+d/2,                   -(1+b)*r + (1-b)*d/2]];

[VSS,DSS] = eig(L); % Eigenvalues and eigenvectors

rho = VSS(:,1)/sum(VSS(:,1));

JSS = L(2,1)*rho(1) - L(1,2)*rho(2);

% Now do the same for the case with the averaged rates
L_av = r*[[-1-b, b, 0, 1];
          [1, -1-b, b, 0];
          [0, 1, -1-b, b];
          [b, 0, 1, -1-b];];

[VSSav,DSSav] = eig(L_av); % Eigenvalues and eigenvectors

rho_av = VSSav(:,1)/sum(VSSav(:,1));

JSS_av = L_av(2,1)*rho_av(1) - L_av(1,2)*rho_av(2);

% Plug in some values and plot
r_val = 1; % arb
b_val = 0.1; % Forward rates twice as fast as reverse rates
d_list = 0.01:0.15:2*r_val;
JSS_plot = subs(JSS,[r,b],[r_val,b_val]);
JSS_plot = double(subs(JSS_plot,d,d_list));

JSS_av_plot = double(subs(JSS_av,[r,b],[r_val,b_val]));

figure; hold on
plot(d_list, JSS_plot)
avline = yline(JSS_av_plot,'--k');
hold off

%% Blocks of length 1 (fast alternation)


% Symbolic manipulation with the rate matrix
syms L1
syms r a b

% r: the average forwards rate between the two regions
% a: ratio of forwards rates between the two regions
% b: the ratio of backwards rate to forwards rate (common between the two regions)

L1 = r*[[-1-a*b, b, 0, a];
       [1, -a-b, a*b, 0];
       [0, a, -1-a*b, b];
       [a*b, 0, 1, -a-b]];

[V1, D1] = eig(L1);

rho1 = V1(:,1)/sum(V1(:,1));

J1 = L1(2,1)*rho1(1) - L1(1,2)*rho1(2);

% Plug in some values and plot
r_val = 1; % arb
b_val = 0.5; % Forward rates twice as fast as reverse rates
a_list = 0.01:0.15:20;

J1_plot = subs(J1,[r,b],[r_val, b_val]);
J1_plot = double(subs(J1_plot,a,a_list));

plot(a_list,J1_plot)
