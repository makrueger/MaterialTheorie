function [sv,sdv] = stress_maxwell(eneu,de,s,sdv,dt,mat)

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

% Elasticity tensor ( see equation 5.14 )
Ee = 2*G*IsymDev + K*IdyI;

Eev = voigt3d(Ee);
dev = voigt3d(de);

% trial stress ( see equation 5.19 )
sigma_trv = voigt3d(s) + Eev*dev;

% see equations 5.20 and 5.21
sv = voigt3d(aG*IsymDev + aK/3*IdyI) * sigma_trv;


end

