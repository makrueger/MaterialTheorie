function [sv,cv,sdvup] = materialbox_vonmises(e,sdv,mat, mflag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3d von Mises plasticity for small strains %
%                                           %
% edit by R.Holtermann 03.07.2014           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% material parameters
emod = mat(1); % Young's modulus
nu   = mat(2); % Poisson ratio
H    = mat(3); % Hardening modulus
sigy = mat(4); % Initial yield stress
r    = mat(5); % r=1: isotropic hardening, r=0: kinematic hardening

% calculate bulk modulus and shear modulus
K = emod/(3*(1-2*nu)); % K, bulk modulus
G = emod/(2*(1+nu));   % G, shear modulus

% internal variable of previous time step
ep = sdv.ep;
k  = sdv.k;
a  = sdv.a;

% unit tensors
I = eye(3);
[I4sym,IdyI] = unitTensor();
I4symdev = I4sym - IdyI/3;

% step 1: calculate trial values
sigvol          = K*(T2_colon_T2(e,I))*I;  %or simpler: sigvol = K*trace(e)*I
edev            = e - trace(e)/3*I;
sigdev_trial    = 2*G*(edev-ep);
kappa_trial     = -r*H*k;
alpha_trial     = -2/3*(1-r)*H*a;
sigdevred_trial = sigdev_trial + alpha_trial;
phi_trial       = norm(sigdevred_trial,'fro') - sqrt(2/3)*(sigy-kappa_trial);

% step 2: check admissible cylinder in deviatoric stress space
if phi_trial <= 0
    sigdev = sigdev_trial;
    sig    = sigvol+sigdev;
    c1     = 2*G;
    c2     = 0;
    n      = zeros(3,3);
    
else
% step 3: radial return
    dlambda = phi_trial/(2*G+2/3*H);
    k       = k+dlambda*sqrt(2/3);
    n       = sigdevred_trial/norm(sigdevred_trial,'fro');
    a       = a+dlambda*n;
    ep      = ep+dlambda*n;
    sigdev  = sigdev_trial- 2*G*dlambda*n;
    sig     = sigvol+sigdev;
    
    c1      = 2*G*( 1-((2*G)*dlambda) / (norm(sigdevred_trial,'fro')) );
    c2      = 4*G^2 *((dlambda)/(norm(sigdevred_trial,'fro'))-(1/(2*G)));   
end

if mflag == 1
    c  = K*IdyI + c1*I4symdev + c2*T2_dyadic_T2(n,n);
else
    c  = moduli_num_vonmises(e,sdv,mat, sig);
end

% return values:
cv = voigt3d(c);
sv = voigt3d(sig);

sdvup.ep = ep;
sdvup.k  = k;
sdvup.a  = a;
end

