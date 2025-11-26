function [u_interior, nodos, time, error] = els_finitos(VhMethod,integrationMethod,N)
a = 0; b = 1;

% f = @(x) (1+pi^2)*sin(pi*x);
% u_exact = @(x) sin(pi*x);

% f = @(x) (1+pi^2)*cos(pi*x)-1;
% u_exact = @(x) cos(pi*x)-1;

u_exact = @(x) (x.^4 - x);
f = @(x) (x.^4 - 12*x.^2 - x);

switch VhMethod
    case 1
        tic
        [u_interior, nodos] = pol_lineales(N,a,b,f,integrationMethod);
        time = toc;
        % Calculo el error fuera para no contaminar tiempo de ejecución
        error = max(abs(u_exact(nodos) - u_interior));
    case 2
        tic
        [u_interior, nodos] = pol_cuadraticos(N,a,b,f,integrationMethod);
        time = toc;
        % Calculo el error fuera para no contaminar el tiempo de ejecución
        % error = max(abs(u_exact(nodos) - u_interior)) ;
        error = sqrt(mean((u_exact(nodos) - u_interior).^2 ));
end
end

function [u_interior, nodos] = pol_lineales(N, a, b, f, integrationMethod)
% Aproximación por metodo de elementos finitos lineales para el
% caso -u''-u = f con nodos equiespaciados

nodos = linspace(a,b,N+1)';
h = nodos(2) - nodos(1);

% Definimos las matrices de rigidez y masa como sparse
% Tienen forma especifica al ser nodos equiespaciados
M = h/6 * spdiags([1,4,1],-1:1,N-1,N-1);
K = 1/h * spdiags([-1,2,-1],-1:1,N-1,N-1);

% El problema es de la forma (M+K)c=F
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

% Solo queremos nodos interiores
F_interior = F(2:end-1);

u_interior = S  \ F_interior;

% Nodos interiores
nodos = nodos(2:end-1);
end

function [u_interior, nodos] = pol_cuadraticos(N, a, b, f, integrationMethod)

phi_1 = @(x) (2*x-1).*(x-1);
phi_2 = @(x) x.*(1-x)*4;
phi_3 = @(x) (2*x-1).*x;

num_nodos = 2*N+1; % Consideramos aquí los nodos intermedios también
nodos = linspace(a,b,num_nodos)';
% El h es el valor entre 2 nodos sin considerar el central, ya que el cambio lo hago desde
% x_k, x_k+1 a 0,1 defino cuadraticos en 0,1 luego el nodo central vendra dado a posteriori
h = nodos(3) - nodos(1); % aprovecho el calculo de la malla

% Las matrices de masa y rigidez (M y K) son pentadiagonales con valores nulos en las (4,2),(2,4)
% hay que hacer el ensamblado teneindo en cuenta el particionado en subdivisión de matrices 3x3
d_main = zeros(num_nodos,1);
d_mid = zeros(num_nodos,1);
d_out = zeros(num_nodos,1);

% Ensamblado de M
d_main(2:2:end) = 16;
d_main(1:2:end) = 8;
d_main(1) = 4; d_main(end) = 4; % Corrección en los bordes

d_mid(1:end) = 2;

d_out(1:2:end) = -1;
d_out(2:2:end) = 0;

M = h/30 * spdiags([d_out d_mid d_main d_mid d_out], -2:2, num_nodos, num_nodos);

% Ensamblado de K
% d_main(2:2:end) = 16; reuso el de M
d_main(1:2:end) = 14;
d_main(1) = 7; d_main(end) = 7;

d_mid(:) = -8;
d_out(1:2:end) = 1; d_out(2:2:end) = 0;

K = (1/(3*h)) * spdiags([d_out d_mid d_main d_mid d_out], -2:2, num_nodos, num_nodos);

S = K + M;

F = zeros(num_nodos,1);

% Iterar sobre el número de elmentos no nodos
for k = 1:N
    % Indexamos relativo al elemento
    id_cent = 2*k;
    id_izq = 2*k-1;
    id_der = 2*k+1;

    % Para el cambio de variable
    x_k = nodos(id_izq);

    f_local_1 = @(x) f(x_k+x*h) .* phi_1(x) * h;
    f_local_2 = @(x) f(x_k+x*h) .* phi_2(x) * h;
    f_local_3 = @(x) f(x_k+x*h) .* phi_3(x) * h;

    int_1 = aprox_integral(integrationMethod,0,1,f_local_1,10);
    int_2 = aprox_integral(integrationMethod,0,1,f_local_2,10);
    int_3 = aprox_integral(integrationMethod,0,1,f_local_3,10);

    % Ensamblar usando los indices reasignados
    F(id_izq) = F(id_izq) + int_1;
    F(id_cent) = F(id_cent) + int_2;
    F(id_der) = F(id_der) + int_3;
end

% Considaramos condición de nodos interiores
S = S(2:end-1,2:end-1);
F = F(2:end-1);

u_interior =  S \ F;
nodos = nodos(2:end-1);

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


