function [u_interior, nodos] = pol_lineales(N, a, b, f, integrationMethod)
nodos = linspace(a,b,N+1)';
h = nodos(2) - nodos(1);

% Definimos las matrices de rigidez y masa como sparse
% Tienen forma especifica al ser nodos fijos
M = h/6 * spdiags([1,4,1],-1:1,N-1,N-1);
K = 1/h * spdiags([-1,2,-1],-1:1,N-1,N-1);

S = M+K;


% Ensamblage de F
F = zeros(N+1,1);
for k = 1:N
    % Hago ensablado de F por cada elemento (x_k, x_{k+1})
    % Considero el cambio de variable x = x_k + e*h que traslada de x en [x_k, x_{k+1}] hacia e en [0,1]


    x_k = nodos(k);

    % Defino handle izquierdo y derecho del producto de f (y los elementos vienen considerando un cambio
    % de intervalo del intervalo de partida al 0,1.

    % Aunque esto sea mas logico según yo resulta que matlab puede considerar handles dentro de handles
    % y elementos fijos (vease x_k) dentro de ellos también
    % f_1 = f(x_k + x * nodos) .* (1 - nodos) * h;
    % f_2 = f(x_k + x * nodos) .* nodos * h;

    % Uso esta versión porque es más facil de implementar los métodos de integración numérica a posteriori
    f_local_1 = @(x) f(x_k + x*h) .* (1 - x) * h;
    f_local_2 = @(x) f(x_k + x*h) .* x * h;

    % Calculo integral local (respecto al elemento k-esimo)

    % Proyección sobre varphi_k (la que baja) en el intervalo [x_k, x_{k+1}]
    int_1 = aprox_integral(integrationMethod,0,1,f_local_1,10);

    % Proyección sobre varphi_k+1 (la que sube) en el intervalo [x_k, x_{k+1}]
    int_2 = aprox_integral(integrationMethod,0,1,f_local_2,10);

    % F(k) = F((k-1)+1) será todos los cachos que involucren en a varphi_k eso son la subida de [x_{k-1}, x_k] (de la iteración previa)
    % Y el aporte de la bajada en [x_k, x_{k+1}]
    F(k) = F(k) + int_1;
    F(k+1) = F(k+1) + int_2;
end

% Solo queremos nodos interiores (los extremos se anulan)
F_interior = F(2:end-1);

u_interior = S  \ F_interior;
end


function [aproxIntegral] = aprox_integral(integrationMethod,a,b,f,n)
% Malla equiespaciada
h = (b-a) / n;
malla = linspace(a,b,n+1); % N subintervalos, N+1 nodos

% Definimos índices para no repetirlos todo el rato
i = 1 : n;       % Lado izquierdo (1:end)
j = 2 : n + 1;   % Lado derecho (2:end)

switch integrationMethod
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
