% L_explicit.m

% A stochastic rate matrix for explicit simiulations with many sites, using
% the same model as for the inifinite random walk. As I have written this
% function, the boundaries of the long chain are reflective - alternatively
% they may be treated as absorbing states, or may include "leaks" leading
% to a gradual decrease in total population.

% Matthew Gerry, April 2023

function [L, sites, block_types] = L_explicit(mA, mB, bias, ga_av, dga, tau, numsites)
    
    if numsites<=mA+mB
        warning("chain length is shorter than one period of the rate variations.")
    end % warning

    k = tau^2./[ga_av + 0.5*dga, ga_av - 0.5*dga]; % Forward rates based on differing gamma
    if bias==Inf
        k_r = [0,0]; % Suppress reverse transitions in infinite bias limit
    else
        k_r = k*exp(-bias); % Reverse rates
    end
    
    sites = (1:numsites) - 0.5*(numsites+1); % Just an array of the site labels, from -(numsites-1)/2 to (numsites-1)/2
    block_types = 1 + ( rem(rem(sites,mA+mB)+mA+mB,mA+mB)>=mA); % 1 for sites in block type A, 2 if in B
    
    L = zeros(numsites); % Pre-allocate rate matrix
    
    for ii=2:numsites-1 % Leave out edge cases for now
        % Populate the matrix with rates - rate type for the transtion
        % determined by site on the left side of the pair
        L(ii,ii+1) = k_r(block_types(ii)); % Reverse rates into site ii from the right
        L(ii,ii-1) = k(block_types(ii-1)); % Forward rates into site ii from the left
        L(ii,ii) = -k(block_types(ii)) - k_r(block_types(ii-1)); % Rates out of site ii
    end % ii

    % Fill in remaining matrix elements
    L(1,2) = k_r(block_types(1));
    L(2,1) = k(block_types(1));

    L(numsites, numsites-1) = k(block_types(numsites-1));
    L(numsites-1, numsites) = k_r(block_types(numsites-1));

    % Populate the last two diagonals
    L(1,1) = -k(block_types(1));
    L(numsites, numsites) = -sum(L(:, numsites));

end % function