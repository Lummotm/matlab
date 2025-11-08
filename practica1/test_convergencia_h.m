% Test de convergencia para verificar convergencia de h
for i = 1:3 
    choiceMethod = i; 
    choiceError = 1; 
    J_values = [200:200:1000]; 
    if i == 1
        N_values = [];
    else
        N_values = J_values; 
    end
    h_values = 1./J_values;
    fprintf("Ejecutando para método %i\n",i)

    % Ejecutamos la función
    [Times, Errors, ~, ~, ~] = practica1_1(choiceMethod, J_values, N_values, choiceError);

    % Extraemos SOLO la diagonal, que tiene mu=0.4 constante
    Errors_diag = diag(Errors);

    % Hacemos polyfit entre log(h) y log(Error) de la diagonal
    % p(1) será la pendiente (el orden p)
    % p(2) será la ordenada en el origen (log(C))
    p = polyfit(log(h_values), log(Errors_diag), 1);

    orden_convergencia = p(1);
    fprintf('El orden de convergencia en espacio (p) es: %.4f\n', orden_convergencia);
    clear all
    close all
end
