% ------------------------------------------------------
% Comparaison NUMÉRIQUE tau_py (Python) vs tau_m (MATLAB)
% ------------------------------------------------------

% 1) Liste des variables symboliques présentes
vars = symvar([tau_py(:); tau_m(:)]);   % toutes les variables qui apparaissent
nvars = numel(vars);

fprintf('Nombre de variables libres : %d\n', nvars);

for trial = 1:5
    % 2) Générer des valeurs aléatoires pour chaque variable
    vals = num2cell(randn(1, nvars));   % gaussien, ou rand pour [0,1]

    % 3) Évaluer les deux versions
    v_py = double(subs(tau_py, vars, vals));
    v_m  = double(subs(tau_m,  vars, vals));

    % 4) Différence max
    diff_max = max(abs(v_py(:) - v_m(:)));
    fprintf('Essai %d : diff max = %.3e\n', trial, diff_max);
end