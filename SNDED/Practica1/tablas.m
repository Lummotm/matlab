%% Generación de Tablas (P1 y P2)
clear; clc;

% Configuración
J_vals_espacio = [50, 100, 200, 400, 800];
J_fijo_tiempo = 2000;
T = 0.5;

%% --- PARTE 1: ECUACIÓN HOMOGÉNEA ---
disp('=== PARTE 1: HOMOGÉNEA ===');

% Tabla 1: Convergencia Espacial (Crank-Nicolson)
fprintf('\n---> Tabla 1: Espacial (CN, mu=0.4)\n');
[~, Err_Mat, ~, ~, ~] = practica1_1(3, J_vals_espacio, []);
imprimir_tabla_espacial(J_vals_espacio, 1./J_vals_espacio, diag(Err_Mat));

% Tabla 2: Convergencia Temporal
fprintf('\n---> Tabla 2: Temporal (J=%d)\n', J_fijo_tiempo);
N_vals_t = [100, 200, 400, 800, 1600];
k_vec = T ./ N_vals_t;

[~, Err_Imp, ~, ~, ~] = practica1_1(2, J_fijo_tiempo, N_vals_t);
[~, Err_CN, ~, ~, ~]  = practica1_1(3, J_fijo_tiempo, N_vals_t);

imprimir_tabla_temporal(N_vals_t, k_vec, Err_Imp, Err_CN);


%% --- PARTE 2: ECUACIÓN NO HOMOGÉNEA ---
disp(' ');
disp('=== PARTE 2: NO HOMOGÉNEA ===');

% Tabla 3: Convergencia Espacial P2 (Crank-Nicolson)
fprintf('\n---> Tabla 3: Espacial P2 (CN, mu=0.4)\n');
[~, Err_Mat_P2, ~, ~, ~] = practica1_2(3, J_vals_espacio, [], T);
imprimir_tabla_espacial(J_vals_espacio, 1./J_vals_espacio, diag(Err_Mat_P2));

% Tabla 4: Convergencia Temporal P2
fprintf('\n---> Tabla 4: Temporal P2 (J=%d)\n', J_fijo_tiempo);
[~, Err_Imp_P2, ~, ~, ~] = practica1_2(2, J_fijo_tiempo, N_vals_t, T);
[~, Err_CN_P2, ~, ~, ~]  = practica1_2(3, J_fijo_tiempo, N_vals_t, T);

imprimir_tabla_temporal(N_vals_t, k_vec, Err_Imp_P2, Err_CN_P2);


%% --- FUNCIONES DE FORMATO ---

function imprimir_tabla_espacial(J, h, Error)
[J, h, Error] = deal(J(:), h(:), Error(:)); % Forzar columnas

Ratio = [NaN; Error(1:end-1) ./ Error(2:end)];
h_ratio = h(1:end-1) ./ h(2:end);
Orden_p = [NaN; log(Ratio(2:end)) ./ log(h_ratio)];

fprintf('%-6s %-10s %-12s %-8s %-8s\n', 'J', 'h', 'Error', 'Ratio', 'Orden p');
fprintf('--------------------------------------------------\n');
for i = 1:length(J)
    fprintf('%-6d %-10.5f %-12.2e %-8.2f %-8.4f\n', J(i), h(i), Error(i), Ratio(i), Orden_p(i));
end
end

function imprimir_tabla_temporal(N, k, Err_Imp, Err_CN)
[N, k, Err_Imp, Err_CN] = deal(N(:), k(:), Err_Imp(:), Err_CN(:)); % Forzar columnas

% Cálculo orden q (Implícito y CN)
k_ratio = k(1:end-1) ./ k(2:end);

Ratio_Imp = [NaN; Err_Imp(1:end-1) ./ Err_Imp(2:end)];
Ord_Imp = [NaN; log(Ratio_Imp(2:end)) ./ log(k_ratio)];

Ratio_CN = [NaN; Err_CN(1:end-1) ./ Err_CN(2:end)];
Ord_CN = [NaN; log(Ratio_CN(2:end)) ./ log(k_ratio)];

fprintf('%-6s %-10s | %-12s %-8s | %-12s %-8s\n', 'N', 'k', 'Err Imp', 'Ord q', 'Err CN', 'Ord q');
fprintf('------------------------------------------------------------------\n');
for i = 1:length(N)
    fprintf('%-6d %-10.5f | %-12.2e %-8.4f | %-12.2e %-8.4f\n', ...
        N(i), k(i), Err_Imp(i), Ord_Imp(i), Err_CN(i), Ord_CN(i));
end
end
