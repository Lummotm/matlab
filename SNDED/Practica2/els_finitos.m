function [] = els_finitos(VhMethod,aproxMethod)
a = 0; b = 1;
u_exact = @(x) sin(pi*x);
f_handle = @(x) (1+pi^2)*sin(pi*x);

% Escoger el espacio Vh (de elementos finitos)
% Mallas uniformes
% La malla depende de la discretización
switch VhMethod
    case 1 % Polinomios lineales
    case 2 % Polinomios cuadráticos
end

end

function [] = pol_lineales(aproxMethod,f,borde_izquierdo,borde_derecho)
N = 100; % Número de elementos finitos

% Aproximación lineal, la derivada es directa, si baja -1 y si sube 1 como estan ordenados pues todo feten
h = nodos(2) - nodos(1); % Aprovechamos la discretización de la malla

% La aproximación no requiere el calculo de las integrales, considerando nodos genericos equiesapciados los productos son fijos
% Se tiene que (u(0)=u(1)= 0) asi que la dimesión de las matrices es N-1
M = h * spdiags([1/6,2/3,1/6],N-1); % Matriz de masa
K = 1/h * spdiags([-1,2,-1],N_1); % Matriz de rigidez



end


function [aproxIntegral] = aprox_integral(aproxMethod,a,b,f,n)
% Malla equiespaciada
h = (b-a) / n;
malla = linspace(a,b,n+1); % N subintervalos, N+1 nodos

% Definimos índices para no repetirlos todo el rato
i = 1 : n;       % Lado izquierdo (1:end)
j = 2 : n + 1;   % Lado derecho (2:end)

switch method
    case 1 % Rectángulo (Izquierda)
        aproxIntegral = h * sum(f(malla(i)));

    case 2 % Punto Medio
        aproxIntegral = h * sum(f((malla(i) + malla(j)) / 2));

    case 3 % Trapecio
        aproxIntegral = (h / 2) * sum(f(malla(i)) + f(malla(j)));

    case 4 % Simpson
        c = (malla(i) + malla(j)) / 2; % Ptos. centrados
        aproxIntegral = (h / 6) * sum(f(malla(i)) + 4*f(c) + f(malla(j)));
end

end
