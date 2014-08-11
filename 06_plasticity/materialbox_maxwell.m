function [s,c,sdv] = materialbox_maxwell(eneu,de,s,sdv,dt,mat,mflag)

% calculation of stresses
[s,sdv] = stress_maxwell(eneu,de,s,sdv,dt,mat);

% calculation of moduli
if mflag == 1
    % analytically
    c = moduli_maxwell(sdv,dt,mat);
elseif mflag == 2
    % numerically
    % to be implemented
    c = moduli_maxwell(sdv,dt,mat);
end;


end

