function [s,sdv] = stress_hooke(ev,sdv,mat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hooke's law for isotropic elasticity %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% material parameters
lambda = mat(1);
mu     = mat(2);

% 2nd order unit tensor
I = eye(3);

% 2nd order unit tensor in Voigt notation
deltav = voigt3d(I);

% cauchy stress
s = lambda*(ev(1) + ev(2) + ev(3))*deltav + 2*mu*ev;

end

