% modular_RW_bias.m

% Statistcs for a classical random walk with transition rates based on a
% quantum mechanical model with tunnelling coupling between sites on a
% molecule/chain and local decoherence. Decoherence rates alternate in
% modules of sizes nA and nB (specified). Here we investigate the case with
% bias (constant ratio of forward to reverse rates passed to the function
% diffusionLchi which builds the counting field-dressed Liouvillian).

% Matthew Gerry, February 2023

% Parameters determining rates
tau = 1; % tunnel coupling between nearest neighbours (all sites degenerate)
ga_av = 10; % average decoherence rate between the two regions
bias = 1.5; % bias factor (ratio of forward to reverse rates)

dchi = 0.01; % Counting field step


%%% PLOT CGF VS CHI FOR ONE CHOICE OF BLOCK LENGTH AND dga

% Parameters specific to CGF-plotting section of the program (arb.)
nA_cgf = 2;
nB_cgf = 3;
dga_cgf = 10;
chisteps_cgf = 61; % Must be odd

% Calculate CGF using diffusionLchi and CGFclassical functions
[Lchi, ~, ~, chi] = diffusionLchi(nA_cgf,nB_cgf,bias,ga_av,dga_cgf,tau,dchi,chisteps_cgf);
CGF_plot = CGFclassical(Lchi);

% Plotting imaginary and real parts of the cumulant generating functions
figure;
subplot(1,2,1); hold on; box on
plot(chi, imag(CGF_plot), DisplayName="$\mathcal{G}(\chi)$")
xlabel("$\chi$",Interpreter="latex")
legend(Interpreter="latex")
set(gca, fontsize=14)
hold off

subplot(1,2,2); hold on; box on
plot(chi, real(CGF_plot), DisplayName="$\mathcal{G}(\chi)$")
xlabel("$\chi$",Interpreter="latex")
legend(Interpreter="latex")
set(gca, fontsize=14)
hold off


%%% FCS FOR A RANGE OF BLOCK LENGTHS AND dga VALUES

dga_list = 0:0.05:1.9*ga_av; % difference between decoherence rates will vary

% Hold nA (size of blocks associated with larger gamma) fixed and vary nB
nA = 3;
nB_list = [1,3,5,11];

% Pre-allocation of lists to contain data - weighted average
J_wa = zeros(length(nB_list),length(dga_list)); % Flux based on "weighted average" rate
                                                %   Now changes with each dga value
S_wa = zeros(length(nB_list),length(dga_list));
C3_wa = zeros(length(nB_list),length(dga_list));                         
C4_wa = zeros(length(nB_list),length(dga_list));

% Pre-allocation of lists to contain data - modular RW
J = zeros(length(nB_list),length(dga_list));
S = zeros(length(nB_list),length(dga_list));
C3 = zeros(length(nB_list),length(dga_list));                         
C4 = zeros(length(nB_list),length(dga_list));

for jj=1:length(nB_list) % Iterate through nB values
    nB = nB_list(jj);
    n = [nA,nB];
    for ii=1:length(dga_list) % Iterate through delta gamma values
        dga = dga_list(ii);

        ga_wa = ga_av - 0.5*dga*diff(n)/sum(n); % Weighted average gamma based on relative block lengths

        % Use diffusionLchi function to derive chi-dressed Liouvillians, rates
        [Lchi, k, kr, ~] = diffusionLchi(nA,nB,bias,ga_av,dga,tau,dchi,5); % modular
        [Lchi_wa, k_wa, kr_wa, chi] = diffusionLchi(nA,nB,bias,ga_wa,0,tau,dchi,5); % homogeneous

        % Full counting statistics
        CGF = CGFclassical(Lchi); % Obtain CGF by passing Lchi to CGFclassical function
        CGF_wa = CGFclassical(Lchi_wa);

        % Obtain the cumulants
        J(jj,ii) = differentiateCGF(CGF,dchi,1);
        S(jj,ii) = differentiateCGF(CGF,dchi,2);
        C3(jj,ii) = differentiateCGF(CGF,dchi,3);
        C4(jj,ii) = differentiateCGF(CGF,dchi,4);

        J_wa(jj,ii) = differentiateCGF(CGF_wa,dchi,1);
        S_wa(jj,ii) = differentiateCGF(CGF_wa,dchi,2);
        C3_wa(jj,ii) = differentiateCGF(CGF_wa,dchi,3);
        C4_wa(jj,ii) = differentiateCGF(CGF_wa,dchi,4);

    end % ii
end % jj


%%% Plotting the cumulants

% First cumulant
figure
for jj=1:length(nB_list)
    subplot(1,length(nB_list),jj); hold on; box on
    plot(dga_list, J(jj,:), '-b', DisplayName="Modular")
    plot(dga_list, J_wa(jj,:), '--r', DisplayName="Homogeneous")
    xlabel("$\Delta\gamma$",Interpreter="latex")
    if jj==1
        ylabel("$\langle J\rangle$",Interpreter="latex")
        legend('Location','north',Orientation='horizontal')
    elseif jj==4
        title(strcat("$n_A\;=\;$",num2str(nA),"$,\;\bar{\gamma}\;=\;$",num2str(ga_av)),Interpreter="latex")
    end
    yl = ylim;
    xl = xlim;
    text(xl(1) + 0.05*diff(xl), yl(1) + 0.5*diff(yl), strcat("$n_B\;=\;$",num2str(nB_list(jj))), Interpreter="latex",FontSize=14)
    set(gca,fontsize=14)
end % jj

% Second cumulant
figure
for jj=1:length(nB_list)
    subplot(1,length(nB_list),jj); hold on; box on
    plot(dga_list, S(jj,:), '-b', DisplayName="Modular")
    plot(dga_list, S_wa(jj,:), '--r', DisplayName="Homogeneous")
    xlabel("$\Delta\gamma$",Interpreter="latex")
    if jj==1
        ylabel("$\langle\langle J^2\rangle\rangle$",Interpreter="latex")
        legend('Location','north',Orientation='horizontal')
    elseif jj==4
        title(strcat("$n_A\;=\;$",num2str(nA),"$,\;\bar{\gamma}\;=\;$",num2str(ga_av)),Interpreter="latex")
    end
    yl = ylim;
    xl = xlim;
    text(xl(1) + 0.05*diff(xl), yl(1) + 0.5*diff(yl), strcat("$n_B\;=\;$",num2str(nB_list(jj))), Interpreter="latex",FontSize=14)
    set(gca,fontsize=14)
end % jj

% Third cumulant
figure
for jj=1:length(nB_list)
    subplot(1,length(nB_list),jj); hold on; box on
    plot(dga_list, C3(jj,:), '-b', DisplayName="Modular")
    plot(dga_list, C3_wa(jj,:), '--r', DisplayName="Homogeneous")
    xlabel("$\Delta\gamma$",Interpreter="latex")
    if jj==1
        ylabel("$\langle\langle J^3\rangle\rangle$",Interpreter="latex")
        legend('Location','north',Orientation='horizontal')
    elseif jj==4
        title(strcat("$n_A\;=\;$",num2str(nA),"$,\;\bar{\gamma}\;=\;$",num2str(ga_av)),Interpreter="latex")
    end
    yl = ylim;
    xl = xlim;
    text(xl(1) + 0.05*diff(xl), yl(1) + 0.5*diff(yl), strcat("$n_B\;=\;$",num2str(nB_list(jj))), Interpreter="latex",FontSize=14)
    set(gca,fontsize=14)
end % jj

% Fourth cumulant
figure
for jj=1:length(nB_list)
    subplot(1,length(nB_list),jj); hold on; box on
    plot(dga_list, C4(jj,:), '-b', DisplayName="Modular")
    plot(dga_list, C4_wa(jj,:), '--r', DisplayName="Homogeneous")
    xlabel("$\Delta\gamma$",Interpreter="latex")
    if jj==1
        ylabel("$\langle\langle J^4\rangle\rangle$",Interpreter="latex")
        legend('Location','north',Orientation='horizontal')
    elseif jj==4
        title(strcat("$n_A\;=\;$",num2str(nA),"$,\;\bar{\gamma}\;=\;$",num2str(ga_av)),Interpreter="latex")
    end
    yl = ylim;
    xl = xlim;
    text(xl(1) + 0.05*diff(xl), yl(1) + 0.5*diff(yl), strcat("$n_B\;=\;$",num2str(nB_list(jj))), Interpreter="latex",FontSize=14)
    set(gca,fontsize=14)
end % jj


% ======================================================================= %
%%% Function defintion

function cumulant = differentiateCGF(CGF,dchi,order)
% Give the values of the cumulants at any order up to one less than the length of the CGF array
    if order >= length(CGF)
        error("Order argument must be at least 1 less than the length of the CGF array.")
    end
    differences = diff(CGF, order); % Finite differences at the order specified
    l = length(differences);
    if rem(l,2)~=0
        diff0 = differences((l + 1)/2); % Finite differences evaluated near chi=0
    elseif rem(l,2)==0
        diff0 = (differences(l/2) + differences(1 + l/2))/2; % Finite differences evaluated near chi=0
        % Averaged values around zero for even-length array
    end
    cumulant = (-1i/dchi)^order*diff0; % Values of a cumulant or specified order
end % function
