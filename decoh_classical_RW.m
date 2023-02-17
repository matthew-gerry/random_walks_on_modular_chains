% decoh_classical_RW.m
% Statistcs for a classical random walk with transition rates based on a
% quantum mechanical model with tunnelling coupling between sites on a
% molecule/chain and local decoherence. Symmetric case - alternating
% decoherence rates but no bias.

% Matthew Gerry, February 2023

% Parameters determining rates
tau = 1; % tunnel coupling between nearest neighbours (all sites degenerate)
ga_av = 10; % average decoherence rate between the two regions

% Counting field step
dchi = 0.005;


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


%%% DO FCS FOR A RANGE OF BLOCK LENGTHS AND dga VALUES

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
