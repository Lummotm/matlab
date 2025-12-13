% tablas_final.m
% Genera tablas de convergencia limpias (EOC) usando 'evalc' para 
% silenciar el output interno de las funciones practica1_1 y practica1_2.

close all; clear; clc;
fprintf('GENERANDO TABLAS DE CONVERGENCIA (SILENCIANDO OUTPUT INTERNO)...\n\n');

% --- CONFIGURACIÓN ---
J_set = [50, 100, 200, 400, 800]; 
N_set = [100, 200, 400, 800, 1600];

% =========================================================================
% TABLA 1: P1 CONVERGENCIA ESPACIAL (h)
% =========================================================================
fprintf('=== TABLA 1: P1 CONVERGENCIA ESPACIAL (h) ===\n');
fprintf('Metodo: Crank-Nicolson\n');
fprintf('%-6s %-12s %-12s %-10s %-10s\n', 'J', 'h', 'Error', 'Ratio', 'Orden(p)');
fprintf('------------------------------------------------------------\n');

h_vals = 1./J_set;
Errors = zeros(length(J_set), 1);

for i = 1:length(J_set)
    J = J_set(i);
    
    % TRUCO: Usamos evalc para ejecutar la orden y "tragarnos" el texto que imprime
    cmd = sprintf('[~, err_mat, ~, ~, ~] = practica1_1(3, %d, []);', J);
    [~] = evalc(cmd); % El output se guarda en ~ (basura), pero 'err_mat' se crea en el workspace
    
    Errors(i) = max(err_mat);
    
    if i > 1
        ratio = Errors(i-1) / Errors(i);
        h_ratio = h_vals(i-1) / h_vals(i);
        p_local = log(ratio) / log(h_ratio);
        fprintf('%-6d %-12.5f %-12.2e %-10.2f %-10.4f\n', ...
            J, h_vals(i), Errors(i), ratio, p_local);
    else
        fprintf('%-6d %-12.5f %-12.2e %-10s %-10s\n', ...
            J, h_vals(i), Errors(i), '-', '-');
    end
end
fprintf('\n');


% =========================================================================
% TABLA 2: P1 CONVERGENCIA TEMPORAL (k) - ERROR ABSOLUTO
% =========================================================================
fprintf('=== TABLA 2: P1 CONVERGENCIA TEMPORAL (k) - Error ABSOLUTO ===\n');
fprintf('%-6s %-10s | %-22s | %-22s\n', '', '', '     IMPLICITO', '   CRANK-NICOLSON');
fprintf('%-6s %-10s | %-10s %-10s | %-10s %-10s\n', 'N', 'k', 'Err(Abs)', 'Ord(q)', 'Err(Abs)', 'Ord(q)');
fprintf('-----------------------------------------------------------------------\n');

J_fixed = 2000;
T = 0.5;
k_vals = T ./ N_set;
Errors_Imp = zeros(length(N_set), 1);
Errors_CN  = zeros(length(N_set), 1);
u_exact_p1 = @(x,t) sin(2*pi*x) .* exp(-4*pi^2 * t);

for i = 1:length(N_set)
    N = N_set(i);
    
    % --- IMPLICITO (Silenciado) ---
    cmd_imp = sprintf('[~, ~, U_int, t_mesh, x_mesh] = practica1_1(2, %d, %d);', J_fixed, N);
    evalc(cmd_imp); 
    
    % Calculo error absoluto manual
    [X, Time] = meshgrid(x_mesh(2:end-1), t_mesh);
    U_ex = u_exact_p1(X', Time');
    Errors_Imp(i) = max(max(abs(U_int - U_ex)));
    
    % --- CRANK-NICOLSON (Silenciado) ---
    cmd_cn = sprintf('[~, ~, U_int, t_mesh, x_mesh] = practica1_1(3, %d, %d);', J_fixed, N);
    evalc(cmd_cn);
    
    % Reusamos U_ex (misma malla)
    Errors_CN(i) = max(max(abs(U_int - U_ex)));
    
    % Impresion
    if i > 1
        k_ratio = k_vals(i-1) / k_vals(i);
        
        ratio_imp = Errors_Imp(i-1) / Errors_Imp(i);
        q_imp = log(ratio_imp) / log(k_ratio);
        
        ratio_cn = Errors_CN(i-1) / Errors_CN(i);
        q_cn = log(ratio_cn) / log(k_ratio);
        
        fprintf('%-6d %-10.5f | %-10.2e %-10.4f | %-10.2e %-10.4f\n', ...
            N, k_vals(i), Errors_Imp(i), q_imp, Errors_CN(i), q_cn);
    else
        fprintf('%-6d %-10.5f | %-10.2e %-10s | %-10.2e %-10s\n', ...
            N, k_vals(i), Errors_Imp(i), '-', Errors_CN(i), '-');
    end
end
fprintf('\n');


% =========================================================================
% TABLA 3: P2 CONVERGENCIA TEMPORAL (k) - ERROR RELATIVO
% =========================================================================
fprintf('=== TABLA 3: P2 CONVERGENCIA TEMPORAL (k) ===\n');
fprintf('%-6s %-10s | %-22s | %-22s\n', '', '', '     IMPLICITO', '   CRANK-NICOLSON');
fprintf('%-6s %-10s | %-10s %-10s | %-10s %-10s\n', 'N', 'k', 'Err(Rel)', 'Ord(q)', 'Err(Rel)', 'Ord(q)');
fprintf('-----------------------------------------------------------------------\n');

T_p2 = 0.5;
k_vals = T_p2 ./ N_set;
Errors_Imp = zeros(length(N_set), 1);
Errors_CN  = zeros(length(N_set), 1);

for i = 1:length(N_set)
    N = N_set(i);
    
    % --- IMPLICITO P2 (Silenciado) ---
    cmd_imp = sprintf('[~, err_mat, ~, ~, ~] = practica1_2(2, %d, %d, %f);', J_fixed, N, T_p2);
    evalc(cmd_imp);
    Errors_Imp(i) = max(err_mat);
    
    % --- CN P2 (Silenciado) ---
    cmd_cn = sprintf('[~, err_mat, ~, ~, ~] = practica1_2(3, %d, %d, %f);', J_fixed, N, T_p2);
    evalc(cmd_cn);
    Errors_CN(i) = max(err_mat);
    
    if i > 1
        k_ratio = k_vals(i-1) / k_vals(i);
        
        ratio_imp = Errors_Imp(i-1) / Errors_Imp(i);
        q_imp = log(ratio_imp) / log(k_ratio);
        
        ratio_cn = Errors_CN(i-1) / Errors_CN(i);
        q_cn = log(ratio_cn) / log(k_ratio);
        
        fprintf('%-6d %-10.5f | %-10.2e %-10.4f | %-10.2e %-10.4f\n', ...
            N, k_vals(i), Errors_Imp(i), q_imp, Errors_CN(i), q_cn);
    else
        fprintf('%-6d %-10.5f | %-10.2e %-10s | %-10.2e %-10s\n', ...
            N, k_vals(i), Errors_Imp(i), '-', Errors_CN(i), '-');
    end
end
fprintf('-----------------------------------------------------------------------\n');
fprintf('RESULTADOS LIMPIOS. Copia estas tablas al Anexo.\n');
