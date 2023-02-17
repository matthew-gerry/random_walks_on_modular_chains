% diffusionLchi.m
% Function to generate the chi-dressed Liovillian for a random walk on a
% line with alternating blocks of sites with varying transition rates
% based on two different values of 1/k. This random walk has no bias.

% The chi-dressed Liouvillian is that of a corresponding Markov jump
% process on a cyclic network of states, doing counting statistics on
% the number of trips the system has taken around this cycle multiplied by
% the steps-per-trip (number of sites)

% Specifically, the transition rates for each block are given by a
% tunnelling parameter squared, tau^2 divided by a decoherence rate ga
% (based on a quantum->classical model).

% This function also outputs the list of rates (length 2) and the list of
% chi values used in case they are needded for further manipulations.

% Matthew Gerry, February 2023


% Parameters
% nA and nB - the lengths (number of sites) of the two blocks with
%             differing transiton rates. Both must be at least 1.
% ga_av     - the average ga value (proportional to reciprocal of
%             transitionrates) associated with the two blocks
% dga       - the difference in ga between the two blocks
% tau       - sqrt of a coefficient by which all rates are scaled
% dchi      - chi_step to use for numerical full counting statistics
% chisteps  - length of list of chi values - a longer list allows one to
%             obtain higher order cumulants. Must be odd and at least 3.

% Pass the average ga and the difference between ga for the two alternating
% regions as arguments to this function.

% The two different rates alternate in blocks of nA and nB sites,
% respectively.



function [Lchi,k,chi] = diffusionLchi(nA,nB,ga_av,dga,tau,dchi,chisteps)

% Parameter checks - break function execution if invalid
if rem(chisteps,2)==0 || chisteps < 3
    fprintf('chisteps argument must be an odd integer, 3 or greater\n\n');
    return;
elseif nA==0 || nB==0
    fprintf('Trivial case - do not use diffusionLchi function\n\n');
    return;
end

% Define rates based on average ga and ga_difference
gaA = ga_av + 0.5*dga;
gaB = ga_av - 0.5*dga;
k = tau^2./[gaA,gaB]; % List of rates

% List of chi values
chi = -0.5*(chisteps-1)*dchi:dchi:0.5*(chisteps-1)*dchi;

dimL = nA + nB; % Dimension of Liouvillian
Lchi = zeros(dimL,dimL,length(chi)); % Initialize Lchi matrix

% Special case - all one-site blocks
if nA==1 && nB==1
    diag_elements = -sum(k)*ones(chisteps,1);
    Lchi12 = k(1)*exp(-dimL*1i*chi) + k(2);
    Lchi21 = k(1)*exp(dimL*1i*chi) + k(2);
    Lchi(1,1,:) = diag_elements;
    Lchi(2,2,:) = diag_elements;
    Lchi(1,2,:) = Lchi12; Lchi(2,1,:) = Lchi21;

else
    L = zeros(dimL); % Initialize the bare Liouvillian matrix

    L(dimL, 1) = k(2); % Rate out of state 1 to the left, by construction
    L(1, dimL) = k(2); % Rate into state 1 from the left, by construction

    sites = 1:dimL; % Just an array of the site labels, 1, 2, ... dimL
    block_types = 2 - (rem(sites,dimL)~=0 & rem(sites,dimL)<=nA); % 1 for sites in block type A, 2 if in B

    % Assign values immediately above the main diagonal of L
    for ii=1:dimL-1
        L(ii,ii+1) = k(block_types(ii));
    end % ii

    % Assign values immediately below the main diagonal of L
    for ii=2:dimL
        L(ii,ii-1) = k(block_types(ii-1)); % ii minus because rate is determined by block type of site to the left
    end % ii

    % Assign diagonal values of L
    for ii=1:dimL
        L(ii,ii) = -sum(L(:,ii)); % Normalization condition
    end % ii

    % Chi-dressed Liouvillian
    for jj=1:chisteps
        Lchi(:,:,jj) = L; % Prepare 3-index object as just a stack of L matrices
    end % jj

    % Incorporate chi-dependence to 1<->2 transitions, scale by number of
    % steps in a cycle (dimL)
    Lchi(1,2,:) = L(1,2)*exp(-dimL*1i*chi);
    Lchi(2,1,:) = L(2,1)*exp(dimL*1i*chi);
    

end % cases





end