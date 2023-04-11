% CGF_detailed.m

% This function derives the cumulant generating function for a classical
% Markov jump process given the counting field-dressed Liouvillian
% The argument is said Liouvillian, a four-index array where the first two
% dimensions are the dimension of the bare Liouvillian and the last two 
% dimensions are the number of points at which we will have the CGF evaluated
% for each of the two counting fields. E.g. Lchi(:,:,ii,jj) is the chi-dressed
% Liouvillian evaluated at particular numerical values of chiA and chiB.

% Output is a 2d array containing values of the CGF at the corresponding
% values of chi at which the Lchi input is evaluated

% Matthew Gerry, April 2023

function CGF = CGF_detailed(Lchi)

    CGF = zeros(size(Lchi,3),size(Lchi,4)); % Initialize CGF list
    for ii=1:size(Lchi,3)
        for jj=1:size(Lchi,4)
            d = eig(Lchi(:,:,ii,jj)); % Obtain eigenvalues of Lchi
            CGF(ii, jj) = d(real(d)==max(real(d))); % Identify dominant eigenvalue as CGF
        end % jj
    end % ii

end % function
