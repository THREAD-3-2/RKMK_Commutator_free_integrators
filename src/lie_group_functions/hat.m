function A = hat(v)
% It returns the skew symmetric matrix A associated to the vector v, such
% that hat(a)b=axb for all 3-component vectors a and b, with "x" the cross product
%
% :param v: vector with 3 components
%
% :returns: element of the Lie algebra so(3) (skew symmetric 3x3 matrix)

    A = zeros(3);
    A(1,2) = -v(3);
    A(1,3) = v(2);
    A(2,3) = -v(1);

    A = A - A';
    
end