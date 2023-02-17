% decoh_classical_RW.m
% Statistcs for a classical random walk with transition rates based on a
% quantum mechanical model with tunnelling coupling between sites on a
% molecule/chain and local decoherence
% Matthew Gerry, February 2023

% Parameters determining rates
tau = 1; % tunnel coupling between nearest neighbours (all sites degenerate)
ga_av = 10; % average decoherence rate between the two regions
dga_list = 0:4:1.9*ga_av; % difference between decoherence rates will vary
% dga_list = [3];

% Counting field
dchi = 0.005;
chi = [-2*dchi,-dchi,0,dchi,2*dchi];

% % SYMMETRIC CASE (alternating decoherence rates but no bias)
% Analytic results for homogeneous chain with averaged decoherence rate
k_av = tau^2/ga_av; % classical transition rates based on model of Cao (NJP 15, 085010, 2013)
% Mean flux at steady state is trivially zero
S_av = 2*k_av; % Classic result for symmetric random walk


block_lengths = 1:6;

S = zeros(length(block_lengths),length(dga_list)); % Initialize lists for variance results
C4 = zeros(length(block_lengths),length(dga_list)); % Initialize lists for 4th cumulant results

for ii=1:length(dga_list)
    for block_length=block_lengths % For each block length
        dga = dga_list(ii);

        % Use diffusionLchi function to derive chi-dressed Liouvillian,
        % rates, counting field
        [Lchi,k,chi] = diffusionLchi(block_length,block_length,ga_av,dga,tau,dchi,5);

        % FCS: determine dominant eigenvalue of Lchi to obtain the CGF
        for jj=1:length(chi)
            [V,D] = eig(Lchi(:,:,jj));
            d = diag(D);
            CGF(jj) = d(real(d)==max(real(d)));
        end
        
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

S_av
S

C4 % There is an issue with the method for block_length=1 - ignore first row