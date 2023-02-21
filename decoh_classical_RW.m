% decoh_classical_RW.m
% Statistcs for a classical random walk with transition rates based on a
% quantum mechanical model with tunnelling coupling between sites on a
% molecule/chain and local decoherence. Symmetric case - alternating
% decoherence rates but no bias.

% Matthew Gerry, February 2023

%%% BASIC PARAMETERS

% Parameters determining rates
tau = 1; % tunnel coupling between nearest neighbours (all sites degenerate)
ga_av = 10; % average decoherence rate between the two regions

% Counting field step
dchi = 0.01;


%%% ANALYTIC RESULTS FOR HOMOGENEOUS CHAIN WITH AVERAGED DECOHERENCE RATE

k_av = tau^2/ga_av; % classical transition rates based on model of Cao (NJP 15, 085010, 2013)
% Mean flux at steady state is trivially zero
S_av = 2*k_av; % Classic result for symmetric random walk


%%% PLOT CGF VS CHI FOR ONE CHOICE OF BLOCK LENGTH AND dga

block_length_CGF = 3;
dga_cgf = 10;
chisteps_cgf = 61; % Must be odd
[Lchi, k, chi] = diffusionLchi(block_length_CGF,block_length_CGF,ga_av,dga_cgf,tau,dchi,chisteps_cgf);
CGF_plot = CGFclassical(Lchi);

figure; hold on; box on
plot(chi, CGF_plot, DisplayName="$\mathcal{G}(\chi)$")
plot(chi, -0.5*S_av*chi.^2,'--',DisplayName="$-\frac{S_0}{2}\chi^2$")
xlabel("$\chi$",Interpreter="latex")
legend(Interpreter="latex")
set(gca, fontsize=14)
hold off


%%% DO FCS FOR A RANGE OF BLOCK LENGTHS AND dga VALUES, EQUAL-SIZED BLOCKS

dga_list = 0:0.05:1.9*ga_av; % difference between decoherence rates will vary
block_lengths = 1:6;

S = zeros(length(block_lengths),length(dga_list)); % Initialize lists for variance results
C4 = zeros(length(block_lengths),length(dga_list)); % Initialize lists for 4th cumulant results

for ii=1:length(dga_list)
    for block_length=block_lengths % For each block length
        dga = dga_list(ii);

        % Use diffusionLchi function to derive chi-dressed Liouvillian,
        % rates, counting field
        [Lchi,k,chi] = diffusionLchi(block_length,block_length,ga_av,dga,tau,dchi,5);

        CGF = CGFclassical(Lchi); % Obtain CGF by passing Lchi to CGFclassical function
        
        % Differentiate the CGF to obtain cumulants
        if block_length==1
            S(1,ii) = 4*prod(k)/sum(k); % Analytic results for one-block sites (J Stat Phys 1454: 1352-1364)
            % Could also do it numerically
        else
            diff2 = diff(CGF,2);
            S(block_length,ii) = -diff2(2)/(dchi^2);
        end

        diff4 = diff(CGF,4);
        C4(block_length,ii) = diff4(1)/(dchi^4);
    
    end % block_lengths
end % ii

figure; hold on; box on
for jj=1:length(block_lengths)
    plot(dga_list, C4(jj,:), DisplayName=strcat("$l\;=\;$",num2str(block_lengths(jj))))
end
xlabel("$\Delta\gamma$",Interpreter="latex")
ylabel("$\langle\langle J^4\rangle\rangle$",interpreter="latex")
legend(Interpreter="latex", location="northwest")
set(gca, fontSize=14)
hold off


%%% FCS WITH UNEVEN-SEZED BLOCKS, NO BIAS

% Hold nA (size of blocks associated with larger gamma) fixed and vary nB
nA = 3;
nB_list = [1,3,5,11];
n = [2,3]; % Fixed block lengths for now, n:=[nA, nB]

% Pre-allocation of lists to contain data
S_wa = zeros(length(nB_list),length(dga_list)); % Variance based on "weighted average" rate
                                                % Now changes with each dga value
C4_wa = zeros(length(nB_list),length(dga_list));
S_ue = zeros(length(nB_list),length(dga_list)); % "ue" - uneven block sizes
C4_ue = zeros(length(nB_list),length(dga_list));

for jj=1:length(nB_list)
    nB = nB_list(jj);
    n = [nA,nB];
    for ii=1:length(dga_list) % Run through list of delta gamma values
        dga = dga_list(ii);
    
        % Use diffusionLchi function to derive chi-dressed Liouvillian,
        % rates, counting field
        [Lchi_ue,k_ue,~] = diffusionLchi(n(1),n(2),ga_av,dga,tau,dchi,5);
        
        % Analytics for homogeneous RW with weighted-averaged gamma
        ga_wa = ga_av - 0.5*dga*diff(n)/sum(n); % Weighted average gamma
        S_wa(jj,ii) = 2*tau^2/ga_wa; % S based on weighted averaged rate (analytic)
        
        % diffusionLchi function for homogeneous RW with weighted average gamma
        [Lchi_wa,k_wa,chi] = diffusionLchi(n(1),n(2),ga_wa,0,tau,dchi,5);
    
        % Full counting statistics
        CGF_ue = CGFclassical(Lchi_ue); % Obtain CGF by passing Lchi to CGFclassical function
        CGF_wa = CGFclassical(Lchi_wa);
    
        % Second and fourth cumulants for alternating chain
        diff2 = diff(CGF_ue,2);
        S_ue(jj,ii) = -diff2(2)/(dchi^2);
        diff4 = diff(CGF_ue,4);
        C4_ue(jj,ii) = diff4(1)/(dchi^4);
    
        % Fourth cumulant for homogeneous chain (we have the second cumulant analytically)
        diff4_wa = diff(CGF_wa,4);
        C4_wa(jj,ii) = diff4_wa(1)/(dchi^4);
    
    end % ii
end % nB

figure
for jj=1:length(nB_list)
    subplot(length(nB_list),2,2*jj-1); hold on; box on
    plot(dga_list, S_ue(jj,:), 'b',DisplayName="Modular")
    plot(dga_list, S_wa(jj,:), '--r',DisplayName="Homogeneous")
    ylabel("$\langle\langle J^2\rangle\rangle$",Interpreter="latex")
    if jj==1
        legend(Location="north",Orientation="horizontal")
    elseif jj==4
        xlabel("$\Delta\gamma$",Interpreter="latex")
    end
    yl = ylim;
    xl = xlim;
    text(xl(1) + 0.05*diff(xl), yl(1) + 0.5*diff(yl), strcat("$n_B\;=\;$",num2str(nB_list(jj))), Interpreter="latex",FontSize=14)
    set(gca,fontsize=14)
    hold off
    
    subplot(length(nB_list),2,2*jj); hold on; box on
    plot(dga_list, C4_ue(jj,:), 'b')
    plot(dga_list, C4_wa(jj,:), '--r')
    ylabel("$\langle\langle J^4\rangle\rangle$",Interpreter="latex")
    if jj==1
        title(strcat("$n_A\;=\;$",num2str(nA),"$,\;\bar{\gamma}\;=\;$",num2str(ga_av)),Interpreter="latex")
    elseif jj==4
        xlabel("$\Delta\gamma$",Interpreter="latex")
    end
    set(gca,fontsize=14)
    hold off
end % jj
