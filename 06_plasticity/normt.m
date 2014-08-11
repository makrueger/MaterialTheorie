function [output] = normt(input)
%NORMT calculates norm of a second order tensor given in index notation
%   for a 2nd order tensor T holds:
%       norm(T) = sqrt( sum( [T_ij]^2) )
    
    temp = 0;
    for i = 1:3
        for j = 1:3
            temp = temp + (input(i,j))^2; 
        end
    end
    output = sqrt(temp);
end

