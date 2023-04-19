% pdf_direct_modRW.m

% Calculate the probability distribution function directly for a modular
% random walk on a long chain, given an initial state (at site zero, by
% choice). Do this by solving the master equation numerically.

% Matthew Gerry, April 2023

function [PDF, sites, n_av, v_av, D_av] = pdf_direct(mA,mB,bias,ga_av,dga,tau,numsites,dt,tmax)

    if rem(numsites,2)==0
        Error("For symmetry, please choose an odd value for the argument numsites")
    end

    p0 = zeros([numsites,1]); p0((numsites+1)/2) = 1; % Initial state (1 at central site)

    % Calculate the rate matrix for this random walk
    [L, sites, ~] = L_explicit(mA,mB,bias,ga_av,dga,tau,numsites);

    % Solve master equation numerically
    time = 0:dt:tmax;

    PDF = zeros(numsites,length(time)); % Pre-allocate time-series of prob dist
    
    for ii=1:length(time)
        t = time(ii);
        PDF(:,ii) = expm(L*t)*p0;
    end

    % Statistics of n 
    dpdt = L*PDF;
    
    n_av = sum(PDF.*repmat(sites',[1,length(time)]));
    v_av = sum(dpdt.*repmat(sites',[1,length(time)]));
    D_av = 0.5*(sum(dpdt.*repmat((sites.^2)',[1,length(time)])) - 2*n_av.*v_av);

end % function