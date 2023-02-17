% CGFclassical.m

% This function derives the cumulant generating function for a classical
% Markov jump process given the counting field-dressed Liouvillian
% The argument is said Liouvillian, a three-index list where the first two
% dimensions are the dimension of the bare Liouvillian and the third
% dimension is the number of points at which we will have the CGF evaluated
% (i.e. dimension of the output). E.g. Lchi(:,:,ii) is the chi-dressed
% Liouvillian evaluated at a particular numerical value of chi.

% Output is a row vector containing values of the CGF at the corresponding
% values of chi at which the Lchi input is evaluated

% Matthew Gerry, February 2023

function CGF = CGFclassical(Lchi)

    CGF = zeros(1,size(Lchi,3)); % Initialize CGF list
    for ii=1:length(CGF)
        d = eig(Lchi(:,:,ii)); % Obtain eigenvalues of Lchi
        CGF(ii) = d(real(d)==max(real(d))); % Identify dominant eigenvalue as CGF
    end % ii

end % function
