% Este test genera la gráfica de eficiencia con todos los nodos de (error,tiempo) de cada método.

close all;
clear;
disp(' INICIANDO TEST DE EFICIENCIA (P1 COMPLETO) ');

% Configuración

[Times_Exp, Errors_Exp,U_int,t,x] = practica1_1(3, 1000,1000);

% 1. Recuperar dimensiones reales de la matriz de resultados
[~, num_cols] = size(U_int);

% 2. Definir los bordes espaciales
% Necesitamos filas de ceros que tengan el mismo ancho que la matriz (num_cols)
borde_x0 = zeros(1, num_cols); % Frontera superior (x=0)
borde_x1 = zeros(1, num_cols); % Frontera inferior (x=1)

% 3. Construir la matriz completa U
% Concatenación VERTICAL: [Fila x=0; Matriz Interior; Fila x=1]
U_full = [borde_x0; U_int; borde_x1];

% 4. Generar la malla para graficar (Meshgrid)
% t y x vienen de tu función practica1_1.
% x tiene dimensión (J+1) y t tiene dimensión (N+1)
[T_grid, X_grid] = meshgrid(t, x);

% 5. Graficar
figure(20);
mesh(T_grid, X_grid, U_full); % O usa surf(T_grid, X_grid, U_full)
xlabel('Tiempo (t)');
ylabel('Espacio (x)');
zlabel('u(x,t)');
title('Solución de la Ecuación del Calor');
colorbar;
