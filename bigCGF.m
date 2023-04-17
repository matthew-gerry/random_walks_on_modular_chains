% bigCGF.m

% This function generates a 4d array representing the CGF where the
% counting field, chi, varies along the fourth index.
% Varying the first three indices amounts to varying the parameters, where
% the first, second, and third indices represent delta_gamma, b, and the
% segment length m respectively. A- and B-segment lengths are equal.
% Note some arguments passed to bigCGF must come in the form of arrays.

% Matthew Gerry, April 2023

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