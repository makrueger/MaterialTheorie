function [I4sym,IdyI] = unitTensor()
I = eye(3);
I4sym = zeros(3,3,3,3);
IdyI  = zeros(3,3,3,3);
for i = 1:3
    for j = 1:3
        for m = 1:3
            for l = 1:3
                I4sym(i,j,m,l) = 0.5*( I(i,m)*I(j,l) + I(i,l)*I(j,m) );
                IdyI(i,j,m,l)  = I(i,j)*I(m,l);
            end
        end
    end
end


end

