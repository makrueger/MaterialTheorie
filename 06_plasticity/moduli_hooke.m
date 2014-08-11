function c = moduli_hooke(ev,sdv,mat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hooke's law for isotropic elasticity %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% material parameters
lambda = mat(1);
mu     = mat(2);

% elasticity tensor

c = [lambda+2*mu lambda      lambda      0    0  0;
     lambda      lambda+2*mu lambda      0    0  0;
     lambda      lambda      lambda+2*mu 0    0  0;
     0           0           0           mu   0  0;
     0           0           0           0    mu 0;
     0           0           0           0    0  mu];


end

