% Funciones base {x, 1-x, 4x-4x^2}
phi =  {
    @(x) x;
    @(x) 1 - x;
    @(x) 4.*x - 4.*x.^2;
    };

der_phi = {
    @(x) 1 + 0.*x;
    @(x) -1 + 0.*x;
    @(x) 4 - 8.*x;
    };

N = length(phi);

M = zeros(N,N);
K = zeros(N,N);

for i = 1:N
    for j = 1:N
        producto = @(x) phi{i}(x) .* phi{j}(x);
        M(i,j) = integral(producto,0,1);

        if (abs(M(i,j)) <= eps)
            M(i,j) = 0;
        end

        producto_derivada = @(x) der_phi{i}(x) .* der_phi{j}(x);
        K(i,j) = integral(producto_derivada,0,1);

        if (abs(K(i,j)) <= eps)
            K(i,j) = 0;
        end
    end
end

format rational
disp("--- Código Final Corregido ---")
disp("Coefs de M (Matriz de Masa)")
disp(M)
disp("Coefs de K (Matriz de Rigidez)")
disp(K)
