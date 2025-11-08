% Test de convergencia para verificar convergencia de k
for i = 2:3 
    choiceMethod = i; 
    choiceError = 1; 
    J_values = [2000]; 
    N_test = [100, 200, 400, 800];
    k_values = T ./ N_test; % T=0.5
    fprintf("Ejecutando para método %i\n",i)

    % Ejecutamos la función
    [Times, Errors, ~, ~, ~] = practica1_1(choiceMethod, J_values, N_values, choiceError);

    % Extraemos SOLO la diagonal, que tiene mu=0.4 constante
    Errors_diag = diag(Errors);

    % Hacemos polyfit entre log(h) y log(Error) de la diagonal
    % p(1) será la pendiente (el orden p)
    % p(2) será la ordenada en el origen (log(C))
    p_k = polyfit(log(k_values), log(Errors(:, 1)), 1);
    orden_q = p_k(1);
    fprintf('El orden de convergencia temporal (q) es: %.4f\n', orden_q);
    clear all
    close all
end
