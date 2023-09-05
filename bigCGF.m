% bigCGF.m

% This function generates a 4d array representing the CGF where the
% counting field, chi, varies along the fourth index. The values of the 
% cumulants at steady state are obtained by differentiating with respect to 
% this index and taking the central value.

% Varying the first three indices amounts to varying the parameters, where
% the first, second, and third indices represent delta_gamma, b, and the
% segment length m respectively. A- and B-segment lengths are equal.

% Note some arguments passed to bigCGF must come in the form of arrays.

% Matthew Gerry, April 2023


% Arguments
% tau         - sqrt of a coefficient by which all rates are scaled
% ga_av       - the average ga value (proportional to reciprocal of
%               transition rates) associated with the two blocks
% dchi        - chi_step to use for numerical full counting statistics
% chisteps    - length of list of chi values - a longer list allows one to
%               obtain higher order cumulants. Must be odd and at least 3
% dga_axis    - list of values of dga (difference in gamma) at which to
%               compute the CGF. This will serve as the horizontal axis of
%               the final plots.
% b_list      - list of values of the bias at which to compute the CGF
% m_list      - list of values of the segment length at which to compute
%               the CGF


function [bigCGF, chi] = bigCGF(tau, ga_av, dchi, chisteps, dga_axis, b_list, m_list)

% Pre-allocate 4D array that will contain CGF at all param values
bigCGF = zeros([length(dga_axis), length(b_list), length(m_list), chisteps]);

% Loop through all sets parameter values and get the CGF as a function of
% chi at each
for ii=1:length(dga_axis)
    dga = dga_axis(ii);
    for jj=1:length(b_list)
        b = b_list(jj);
        for kk=1:length(m_list)
            m = m_list(kk);

            % Calculate CGF from rate matrix using functions in this project
            [Lchi,~,~,chi] = diffusionLchi(m,m,b,ga_av,dga,tau,dchi,chisteps);
            CGF = CGFclassical(Lchi);

            bigCGF(ii,jj,kk,:) = CGF;

        end % kk
    end % jj
end % ii

end % function