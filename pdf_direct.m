% pdf_direct_modRW.m

% Calculate the probability distribution function directly for a modular
% random walk on a long chain, given an initial state (at site zero, by
% choice). This is achieved by solving the master equation numerically. The
% output is a 2-index object where each row is the probability distribution
% at a particular moment in time. Each row represents the distribution one
% time step after the previous row.

% Matthew Gerry, April 2023


% Arguments
% mA        - segment length for A sites
% mB        - segment length for B sites
% bias      - log-ratio of forward to reverse rates (set to zero for
%             unbiased random walk)
% ga_av     - the average ga value (proportional to reciprocal of
%               transition rates) associated with the two blocks
% dga       - the difference in ga between the two blocks
% tau       - sqrt of a coefficient by which all rates are scaled
% numsites  - total number of sites included in the simulation
% dt        - timestep for the evolution of the probability distribution
% tmax      - final time for evolution

function [PDF, sites, n_av, v_av, D_av, C3, C4] = pdf_direct(mA,mB,bias,ga_av,dga,tau,numsites,dt,tmax)

    % Return an error if an even number of sites is used (for consistency
    % with other commands, functions)
    if rem(numsites,2)==0
        Error("For symmetry, please choose an odd value for the argument numsites")
    end
    
    p0 = zeros([numsites,1]); p0((numsites+1)/2) = 1; % Initial state (1 at central site)

    % Calculate the rate matrix for this random walk
    [L, sites, ~] = L_explicit(mA,mB,bias,ga_av,dga,tau,numsites);

    % Solve master equation numerically
    time = 0:dt:tmax;

    PDF = zeros(numsites,length(time)); % Pre-allocate time-series of prob dist
    
    if bias<=0.5 % At low bias, diagonalize L to calculate PDF faster
        [V,D] = eig(L);
        for ii=1:length(time)
            t = time(ii);
            PDF(:,ii) = real(V*expm(D*t)*(V\p0)); % Exponential of L*t acting on p0
        end
    else % Except don't do this at high bias - leads to undesired numerical effects since L is near singular
        for ii=1:length(time)
            t = time(ii);
            PDF(:,ii) = expm(L*t)*p0; % Exponential of L*t acting on p0
        end
    end

    % Use the probability distribution over time to get the time evolution of the cumulants of n 
    dpdt = L*PDF; % Time derivative of the pdf as a function of time
    
    % Mean
    n_av = sum(PDF.*repmat(sites',[1,length(time)]));
    v_av = sum(dpdt.*repmat(sites',[1,length(time)])); % Mean "velocity"
    
    % Diffusion coefficient
    [n_av_grid, sites_grid] = meshgrid(n_av,sites);

    S = sum(PDF.*(sites_grid-n_av_grid).^2); % Second cumulant
    D_av = 0.5*(sum(dpdt.*repmat((sites.^2)',[1,length(time)])) - 2*n_av.*v_av);

    % Skewness (numerical differentiation)
    skw = sum(PDF.*(sites_grid-n_av_grid).^3);
    C3 = diff(skw,1,2)/dt; % Scaled skewness

    % Kurtosis (numerical differentiation)
    krt = sum(PDF.*(sites_grid-n_av_grid).^4) - 3*S.^2;
    C4 = diff(krt,1,2)/dt; % Scaled kurtosis

end % function