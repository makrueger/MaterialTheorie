function c = moduli_num_hooke_RH(ev,sdv,mat,s)

%-------------------------------------------------------------------------
% returns tangent operator associated with Cauchy-stress
% (numerical (consistent) tangent operator conjugate to Cauchy-stress
% stress, for details: see Miehe, CMAME, 134, 223-240 (1996))
% for small strains
%-------------------------------------------------------------------------

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
    [sp,dummy] = stress_hooke(ep,sdv,mat);
    
    % 3. calculation of moduli (forward differential)
    c(:,i) = (sp(:,1) - s(:,1))/eps;

end

