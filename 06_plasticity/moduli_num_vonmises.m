function c = moduli_num_vonmises(e,sdv,mat,s)

%-------------------------------------------------------------------------
% numerical tangent operator for a 2nd order tensor function by 
% forward differential 
%
% input and output defined in "full" tensor notation
%
% Raphael Holtermann 09.07.2014
%-------------------------------------------------------------------------

% perturbation parameter (set to sqrt(computer precision) ideally)
h = 1e-6;

% preallocate tangent tensor (4th order)
c = zeros(3,3,3,3);


for k = 1:3
    for l = 1:3
        % init perturbed strain tensor 
        ep = e;
        
        % perturbate strain tensor (ep)
        ep(k,l) = ep(k,l) + h;
        
        % calculate perturbed stresses
        sp = stress_vonmises(ep,sdv,mat);

        for i = 1:3
            for j = 1:3
                % calculate tangent by forward differential
                c(i,j,k,l) = ( sp(i,j) - s(i,j) ) / h;
            end
        end
    end
end


end %function

