function [s,sdv] = stress_kauderer(ev,sdv,mat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constitutive law of Kauderer %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% material parameter
K      = mat(1);
G      = mat(2);
kapvol = mat(3);
kapdev = mat(4);

c20 = 1/2*K-1/3*G;
c30 = 1/9*K*kapvol;
c40 = 2/27*G*kapdev;
c01 = G;
c02 = 2/3*G*kapdev;
c21 = -4/9*G*kapdev;

% 2nd order unit tensor
I = eye(3);

% 2nd order unit tensor in Voigt notation
deltav = voigt3d(I);

% basic invariants
I1 = ev(1) + ev(2) + ev(3);
I2 = ev(1)^2 + ev(2)^2 + ev(3)^2 + 2*ev(4)^2 + 2*ev(5)^2 + 2*ev(6)^2;

% cauchy stress
s = (2*c20*I1 + 3*c30*I1^2 + 4*c40*I1^3 + 2*c21*I1*I2)*deltav + ...
    (2*c01 + 4*c02*I2 + 2*c21*I1^2)*ev;

end

