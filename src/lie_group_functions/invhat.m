function v = invhat(A)
% Inverse of the hat function
%
% :param A: element of the Lie algebra so(3) (skew symmetric 3x3 matrix)
%
% :returns: 3x1 column vector v

    v = zeros(3,1);
    v(1) = - A(2,3);
    v(2) = A(1,3);
    v(3) = - A(1,2);

end