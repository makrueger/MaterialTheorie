function c = moduli_kauderer(ev,sdv,mat)

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


% 2nd order unit tensor in Voigt notation
deltav = voigt3d(I);

% 4th order unit tensors in Voigt notation
I4symv = voigt3d(I4sym);
IdyIv  = voigt3d(IdyI);

% basic invariants
I1 = ev(1) + ev(2) + ev(3);
I2 = ev(1)^2 + ev(2)^2 + ev(3)^2 + 2*ev(4)^2 + 2*ev(5)^2 + 2*ev(6)^2;

% elasticity tensor
c = (2*c20 + 6*c30*I1 + 12*c40*I1^2 + 2*c21*I2)*IdyIv +...
    (2*c01 + 4*c02*I2 + 2*c21*I1^2)*I4symv +...
    4*c21*I1*(deltav*ev' + ev*deltav')+...
    8*c02*(ev*ev');

end

