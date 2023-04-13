% firstorderCGF_symb.m

% Given a symbolic chi-dressed Liouvillian, perform the symbolic
% manipulations to get the first-order cumulant generating function, as
% defined in Phys Rev 97, 052154. This is given by the determinant of the
% chi-dressed Liouvillians divided by the sum over its first cofactors
% associated with removal of same-indexed rows and columns.

function [dtmnt, cofsum, G1] = firstorderCGF_symb(Lchi)
    dtmnt = det(Lchi); % Determinant
    
    cofsum = 0; % Sum of relevant cofactors
    for ii=1:size(Lchi,1)
        tempmat = Lchi;
        % Remove i'th row and column from matrix
        tempmat(:,ii)= [] ; tempmat(ii,:) = [];

        cof = det(tempmat);
        cofsum = cofsum + cof; % Add cofactor to sum
    end % ii

    G1 = dtmnt/cofsum; % First-order truncated CGF
end % function