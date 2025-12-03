% Calculo de error L2 global
% Uso la misma lógica que con el la forma de calcular el valor, trabajando elemento a elemento
function error = calc_error_L2(VhMethod,nodos,u_h,u_exact)

switch VhMethod
    case 1
        error = error_lineal(nodos,u_h,u_exact);
    case 2
        error  = error_cuadratico(nodos,u_h,u_exact);
end

end

function error = error_lineal(nodos,u_h,u_exact)
phi = {
    @(x) (1-x);
    @(x) (x);
    };

h = nodos(2) - nodos(1);
N = length(nodos) - 1;
error = 0;

for k = 1:N
    x_k = nodos(k);

    u_izq = u_h(k);
    u_der = u_h(k+1);

    % Uso e pq e esta en [0,1] y x en principio esta en el elemento que sea
    u_local = @(e) (u_izq * phi{1}(e) + u_der * phi{2}(e));

    f_integral = @(e) (u_exact(x_k + e*h) - u_local(e)).^2;

    error_local = integral(f_integral,0,1) * h;

    error = error_local + error;
end

error = sqrt(error); % El error global sera la raiz cuadrada de los errores de cada proyección
end


function error = error_cuadratico(nodos,u_h,u_exact)
% La función que aproxima bien la solución será u_h que no se expresa como la suma de los u_h(i)*phi_i(x)
% Uso la misma lógica, que para el ensamblado considero, trabajando por elementos y calculando la contibución de cada uno
% Y usando esta función calculo su aprox en norma L2 (pasando al intervalo de referencia)

phi = {
    @(x) (2*x-1).*(x-1);
    @(x) x.*(1-x)*4;
    @(x) (2*x-1).*x;
    };

N = (length(nodos) - 1) / 2; % Hay 2N+1 nodos, N elementos
h = nodos(3) - nodos(1);
error = 0;

for k = 1:N
    k_izq = 2*k - 1;
    k_cen = 2*k ;
    k_der = 2*k + 1;

    u_izq = u_h(k_izq); u_cen = u_h(k_cen); u_der = u_h(k_der);

    x_k = nodos(k_izq);

    % Uso e pq e esta en [0,1] y x en principio esta en el elemento que sea
    u_local = @(e) (u_izq * phi{1}(e) + u_cen * phi{2}(e) + u_der * phi{3}(e));

    f_integral = @(e) (u_exact(x_k + e*h) - u_local(e)).^2;

    error_local = integral(f_integral,0,1) * h;

    error = error_local + error;
end

error = sqrt(error); % El error global sera la raiz cuadrada de los errores de cada proyección
end
