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
C4 = zeros(length(block_lengths),length(dga_list));

for ii=1:length(dga_list)
    ga = [ga_av - 0.5*dga_list(ii), ga_av + 0.5*dga_list(ii)]; % Higher and lower decoherence rates
    k = tau^2./ga;

    % Numerics for block_lengths greater than 1
    for block_length=block_lengths
        dimL = 2*block_length; % Dimensionality of the Liouvillian
        
        % Define the Liouvillian matrix for the system
        L = zeros(dimL); % Initialize the Liouvillian matrix
        % Populate the entries of L above the main diagonal
        L(1,dimL) = k(2); % By construction, the first site is always at the beginning of a new type-1 block
        for jj=1:dimL-1
            if rem(floor((jj-1)/block_length),2)==0 % Condition to be a site in the A blocks
                L(jj,jj+1) = k(1);
            elseif rem(floor((jj-1)/block_length),2)==1 % Condition to be a site in the B blocks
                L(jj,jj+1) = k(2);
            end
        end % jj
    
        % Populate the entries of L below the main diagonal
        L(dimL, 1) = k(2); % By construction, the first site is always at the beginning of a new type-1 block
        for jj=2:dimL
            if rem(floor((jj-2)/block_length),2)==0 % Condition for the site to the left to be in the A blocks
                L(jj,jj-1) = k(1);
            elseif rem(floor((jj-2)/block_length),2)==1 % Condition for the site to the left to be in the B blocks
                L(jj,jj-1) = k(2);
            end
        end % jj
    
        for jj=1:dimL % Populate diagonals of Liouvillian
            L(jj,jj) = -sum(L(:,jj));
        end % jj
        
        % Full counting statistics
        Lchi = zeros(dimL,dimL,length(chi)); % Initialize counting dependent Liouvillian
        CGF = zeros(1,length(chi)); % Initialize the CGF
        for jj=1:length(chi)
            Lchi(:,:,jj) = L;
        end
        % Just count net steps taken from site 1 to 2
        % Multiply counting field by dimL to account for the number of steps
        % taken around the cycle
        Lchi(1,2,:) = L(1,2)*exp(-dimL*1i*chi);
        Lchi(2,1,:) = L(2,1)*exp(dimL*1i*chi);
    
        for jj=1:length(chi)
            [V,D] = eig(Lchi(:,:,jj));
            d = diag(D);
            CGF(jj) = d(real(d)==max(real(d)));
        end
        
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