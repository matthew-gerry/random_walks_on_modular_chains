% pdf_direct_modRW.m

% Calculate the probability distribution function directly for a modular
% random walk on a long chain, given an initial state (at site zero, by
% choice). Do this by solving the master equation numerically.

% Matthew Gerry, April 2023

function [PDF, sites, n_av, v_av, D_av, C3, C4] = pdf_direct(mA,mB,bias,ga_av,dga,tau,numsites,dt,tmax)

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
    else % Except don't do this at high bias - leads to strange numerical effects since L is near singular
        for ii=1:length(time)
            t = time(ii);
            PDF(:,ii) = expm(L*t)*p0; % Exponential of L*t acting on p0
        end
    end

    % Statistics of n 
    dpdt = L*PDF;
    
    % Mean
    n_av = sum(PDF.*repmat(sites',[1,length(time)]));
    v_av = sum(dpdt.*repmat(sites',[1,length(time)])); % Mean "velocity"
    
    % Diffusion coefficient
    [n_av_grid, sites_grid] = meshgrid(n_av,sites);

    S = sum(PDF.*(sites_grid-n_av_grid).^2);
    D_av = 0.5*(sum(dpdt.*repmat((sites.^2)',[1,length(time)])) - 2*n_av.*v_av);

    % Skewness   
    skw = sum(PDF.*(sites_grid-n_av_grid).^3);
%     C3 = sum(dpdt.*repmat((sites.^3)',[1,length(time)])) - 2*n_av.*sum(dpdt.*repmat((sites.^2)',[1,length(time)])) - 2*sum(PDF.*sites_grid.^2).*v_av + 3*n_av.^2.*v_av;
    C3 = diff(skw,1,2)/dt; % Scaled skewness

    % Kurtosis
    krt = sum(PDF.*(sites_grid-n_av_grid).^4) - 3*S.^2;
%     C4 = sum(dpdt.*repmat((sites.^4)',[1,length(time)])) - 4*sum(dpdt.*repmat((sites.^3)',[1,length(time)])).*n_av - 6*sum(dpdt.*repmat((sites.^2)',[1,length(time)])).*(sum(PDF.*sites_grid.^2) - n_av.^2) - 4*v_av.*(sum(PDF.*sites_grid.^3) - 6*n_av.*sum(PDF.*sites_grid.^2) + 6*n_av.^3);
    C4 = diff(krt,1,2)/dt; % Scaled kurtosis

end % function