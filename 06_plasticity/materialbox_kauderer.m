function [s,c,sdv] = materialbox_kauderer(e,sdv,mat,mflag)

% convert strain tensor in Voigt-format
ev = voigt3d(e);

% calculation of stresses
[s,sdv] = stress_kauderer(ev,sdv,mat);

% calculation of moduli
if mflag == 1
    % analytically
    c = moduli_kauderer(ev,sdv,mat);
elseif mflag == 2
    % numerically
    c = moduli_num_kauderer(ev,sdv,mat,s);
end;

end

