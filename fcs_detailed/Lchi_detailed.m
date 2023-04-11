% Lchi_detailed.m

% Function to generate the chi-dressed Liovillian for a random walk on a
% line with alternating blocks of sites with varying transition rates
% based on two different values of 1/k. We use two different counting
% fields to distinguish between steps in the A and B regions

% The chi-dressed Liouvillian is that of a corresponding Markov jump
% process on a cyclic network of states, doing counting statistics on
% the number of trips the system has taken around this cycle multiplied by
% the steps-per-trip (number of sites)

% Specifically, the transition rates for each block are given by a
% tunnelling parameter squared, tau^2 divided by a decoherence rate ga
% (based on a quantum->classical model).

% This function also outputs the list of rates (length 2) and the 2d-array of
% chi values used in case they are needed for further manipulations.

% Matthew Gerry, February 2023


% Arguments
% nA and nB   - the lengths (number of sites) of the two blocks with
%               differing transiton rates. Both must be at least 1
% bias        - log ratio of forward to reverse rates (set to 0 if no bias)
% ga_av       - the average ga value (proportional to reciprocal of
%               transition rates) associated with the two blocks
% dga         - the difference in ga between the two blocks
% tau         - sqrt of a coefficient by which all rates are scaled
% dchi        - chi_step to use for numerical full counting statistics
% chisteps    - length of list of chi values for each coutning field - a
%               longer list allows one to obtain higher order cumulants.
%               Must be odd and at least 3

% Pass the average ga and the difference between ga for the two alternating
% regions as arguments to this function.

% The two different rates alternate in blocks of nA and nB sites,
% respectively.


function [Lchi, k, k_r, chi] = Lchi_detailed(nA,nB,bias,ga_av,dga,tau,dchi,chisteps)

    % Parameter checks - break function execution if invalid
    if rem(chisteps,2)==0 || chisteps < 3
        error('chisteps argument must be an odd integer, 3 or greater');
    elseif nA==0 || nB==0
        error('Trivial case - use diffusionLchi function only when nA and nB are both nonzero');
    end
    
    bias_ratio = exp(bias); % Ratio of forward to reverse rates

    % Define rates based on average ga and ga_difference
    gaA = ga_av + 0.5*dga;
    gaB = ga_av - 0.5*dga;
    k = tau^2./[gaA,gaB]; % List of rates
    k_r = k/bias_ratio; % List of reverse rates
    
    % List of chi values
    chi = -0.5*(chisteps-1)*dchi:dchi:0.5*(chisteps-1)*dchi;
    [chiA, chiB] = meshgrid(chi, chi);


    dimL = nA + nB; % Dimension of Liouvillian
    Lchi = zeros(dimL,dimL,length(chi),length(chi)); % Initialize Lchi matrix

    % Special case - all one-site blocks
    if nA==1 && nB==1
        % Enter each matrix element directly, vary along chi-dimension as needed
        Lchi(1,1,:,:) = -(k(1) + k_r(2))*ones(chisteps);
        Lchi(2,2,:,:) = -(k_r(1) + k(2))*ones(chisteps);
        Lchi(1,2,:,:) = k_r(1)*exp(-1i*chiA) + k(2)*exp(1i*chiB);
        Lchi(2,1,:,:) = k(1)*exp(1i*chiA) + k_r(2)*exp(-1i*chiB);
    
    else % nA > 1 or nB > 1, dimL > 2
        L = zeros(dimL); % Initialize the bare Liouvillian matrix
    
        L(dimL, 1) = k_r(2); % Rate out of state 1 to the left (reverse direction), by construction
        L(1, dimL) = k(2); % Rate into state 1 from the left (forward direction), by construction
    
        sites = 1:dimL; % Just an array of the site labels, 1, 2, ... dimL
        block_types = 2 - (rem(sites,dimL)~=0 & rem(sites,dimL)<=nA); % 1 for sites in block type A, 2 if in B
    
        % Assign values immediately above the main diagonal of L
        for ii=1:dimL-1
            L(ii,ii+1) = k_r(block_types(ii)); % Reverse transitions
        end % ii
    
        % Assign values immediately below the main diagonal of L
        for ii=2:dimL
            L(ii,ii-1) = k(block_types(ii-1)); % Forward transitions
            % ii-1 because rate is determined by block type of site to the left (by definition)
        end % ii
    
        % Assign diagonal values of L
        for ii=1:dimL
            L(ii,ii) = -sum(L(:,ii)); % Normalization condition
        end % ii
    
        % Chi-dressed Liouvillian
        Lchi = repmat(L,[1,1,chisteps,chisteps]); % Prepare 4-index object as just 2d a stack of L matrices

        % Incorporate the chi-dependence, counting a step at every transition
        Lchi(1,dimL,:,:) = L(1,dimL)*exp(1i*chiB);
        Lchi(dimL,1,:,:) = L(dimL,1)*exp(-1i*chiB);

        for ii=1:dimL-1
            if rem(ii,dimL) <= nA
                Lchi(ii+1,ii,:,:) = L(ii+1,ii)*exp(1i*chiA);
                Lchi(ii,ii+1,:,:) = L(ii,ii+1)*exp(-1i*chiA);
            else
                Lchi(ii+1,ii,:,:) = L(ii+1,ii)*exp(1i*chiB);
                Lchi(ii,ii+1,:,:) = L(ii,ii+1)*exp(-1i*chiB);
            end % cases
        
    
    end % cases
end % function