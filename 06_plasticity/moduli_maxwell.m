function c = moduli_maxwell(sdv,dt,mat)

%%%%%%%%%%%%%%%%%%%%%%%%
% linear Maxwell model %
%%%%%%%%%%%%%%%%%%%%%%%%

% material parameter ( see programming task )
K     = mat(1);
G     = mat(2);
kappa = mat(3);
mu    = mat(4);

% natural relaxation times ( see equation 5.22 )
tG = mu/G;
tK = kappa/K;
aG = 1/(1 + dt/tG);
aK = 1/(1 + dt/tK);

% 2nd order unit tensor
I = eye(3);

% 4th order unit tensors
I4sym = zeros(3,3,3,3);
IdyI  = zeros(3,3,3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                I4sym(i,j,k,l) = 0.5*( I(i,k)*I(j,l) + I(i,l)*I(j,k) );
                IdyI(i,j,k,l)  = I(i,j)*I(k,l);
            end;
        end;
    end;
end;
IsymDev = I4sym - 1/3*IdyI;

% see equation 5.23
ct = 2*G*aG*IsymDev + K*aK*IdyI;

c  = voigt3d(ct);


end

