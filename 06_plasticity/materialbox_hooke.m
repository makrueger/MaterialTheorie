function [s,c,sdv] = materialbox_hooke(e,sdv,mat,mflag)

% convert strain tensor in Voigt-format
ev = voigt3d(e);

% calculation of stresses
[s,sdv]=stress_hooke(ev,sdv,mat);

% calculation of moduli
if mflag==1
    % analytically
    c = moduli_hooke(ev,sdv,mat);
elseif mflag==2
    % numerically
    c = moduli_num_hooke(ev,sdv,mat,s);
end;


end

