function c = moduli_num_maxwell(eneu,de,s,sdv,dt,mat)

%-------------------------------------------------------------------------
% returns tangent operator associated with Cauchy-stress
% (numerical (consistent) tangent operator conjugate to Cauchy-stress
% stress, for details: see Miehe, CMAME, 134, 223-240 (1996))
% for small strains
%-------------------------------------------------------------------------

ev = voigt3d(eneu);

% perturbation parameter
eps = 1e-6;

% preallocation of c-matrix
c = zeros(6,6);

for i=1:6 % calculate tangent columnwise
    % init perturbed strain tensor (voigt notation)
    ep = ev;
    
    % 1. calculation of perturbed strain tensor (ep)
    ep(i) = ep(i) + eps;
    
    % 2. calculation of pertubed stresses (sp, in Voigt-notation)
    [sp,dummy] = stress_maxwell(ep,de,s,sdv,dt,mat);

    % 3. columnwise calculation of moduli
    c(:,i) = (sp(:,1) - s(:,1))/eps;

end

