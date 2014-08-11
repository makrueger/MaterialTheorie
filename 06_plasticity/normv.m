function [output] = normv(input)
%NORMV calculates norm of a second order tensor given in Voigt notation
%   for a 2nd order tensor T holds:
%       norm(T) = sqrt( sum( [T_ij]^2) )
    
    temp = 0;
    %symmetric elements
    for i = 1:3
       temp = temp + (input(i))^2; 
    end
    %skew symmetric elements (have to be summed up twice)
    for i = 4:6
       temp = temp + 2*(input(i))^2; 
    end
    output = sqrt(temp);
end

