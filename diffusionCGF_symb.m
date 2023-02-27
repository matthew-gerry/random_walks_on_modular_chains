%diffusionCGF_symb.m

% A function to generate the scaled cumulant generating function at steady
% state for a unicyclic network representing an random walk on an infinite
% modular chain (classical). The model is the same as that in the
% diffusionLchi.m function, but here we use Matlab's symbolic toolkit for 
% analytic expressions.

% Matthew Gerry, February 2023

% Arguments
% nA and nB   - the lengths (number of sites) of the two blocks with
%               differing transiton rates. Both must be at least 1
% b           - ratio of forwards to reverse rates (set to 1 if no bias)

function G = diffusionCGF_symb(nA, nB)

    if nA==0 || nB==0
        error('Trivial case - use diffusionCGF_symb function only when nA and nB are both nonzero')
    end
    
    % Define variables to feature in expressions for rates, counting field
    % forward rates given by tau^2/(ga_av +/- dga/2)
    syms tau ga_av dga b chi
    
    k = tau^2./[ga_av + 0.5*dga, ga_av - 0.5*dga]; % Forward rates based on differing gamma
    k_r = k/b; % Reverse rates
    
    
    dimL = nA + nB; % Dimension of the rate matrix (number of sites in the cycle)
    Lchi = sym(zeros(dimL)); % Pre-allocate symbolic matrix
    
    % Special case - all blocks one-site-long
    if nA==1 && nB==1
        % Enter each matrix element directly, vary along chi-dimension as needed
        Lchi(1,1) = -(k(1) + k_r(2));
        Lchi(2,2) = -(k_r(1) + k(2));
        Lchi(1,2) = k_r(1)*exp(-dimL*1i*chi) + k_r(2);
        Lchi(2,1) = k(1)*exp(dimL*1i*chi) + k(2);
    
    else
        % Enter the bare rates first
        Lchi(dimL, 1) = k_r(2); % Rate out of state 1 to the left (reverse direction), by construction
        Lchi(1, dimL) = k(2); % Rate into state 1 from the left (forward direction), by construction
    
        sites = 1:dimL; % Just an array of the site labels, 1, 2, ... dimL
        block_types = 2 - (rem(sites,dimL)~=0 & rem(sites,dimL)<=nA); % 1 for sites in block type A, 2 if in B
        
        % Assign values immediately above the main diagonal of L
        for ii=1:dimL-1
            Lchi(ii,ii+1) = k_r(block_types(ii)); % Reverse transitions
        end % ii
    
        % Assign values immediately below the main diagonal of L
        for ii=2:dimL
            Lchi(ii,ii-1) = k(block_types(ii-1)); % Forward transitions
            % ii-1 because rate is determined by block type of site to the left (by definition)
        end % ii
        
        % Assign diagonal values of L
        for ii=1:dimL
            Lchi(ii,ii) = -sum(Lchi(:,ii)); % Normalization condition
        end % ii
    
        % Add counting factors to Lchi where needed
        Lchi(1,2) = Lchi(1,2)*exp(-dimL*1i*chi);
        Lchi(2,1) = Lchi(2,1)*exp(dimL*1i*chi);
    
    end % cases
    
    % Find dominant eigenvalues of Lchi
    d = eig(Lchi);
    % Identify index of dominant eigenvalue
    CGFindex = find(max(double(subs(d,[tau,ga_av,dga,chi],[1.0,1.0,0.5,0.0])))); % Arb. values subbed for tau, ga_av, dga
    
    % Scaled cumulant generating function at steady state
    G = d(CGFindex);


end % function