I = 10;
% 10 da buenos errores, 100 da errores maximales y es el ultimo que ayuda a medir el orden a partir de alli siempre es 10^-12 luego no tiene sentido
% considerar más ya que la pendiente tiende a 1 al ser siempre 10^-12 (epsilon de MATLAB)
J = I*10;
N = I:I:J;
h = 1./N;

time_1 = zeros(1,length(N));
time_2 = zeros(1,length(N));
error_1 = zeros(1,length(N));
error_2 = zeros(1,length(N));

metodos = { 'Rectángulo', 'Pto. Medio', 'Simpson' };


for j = 1:3

    typeError = j;
    for i = 1:length(N)
        [~,time_1(i),error_1(i)] = els_finitos(1,typeError,N(i));
        [~,time_2(i),error_2(i)] = els_finitos(2,typeError,N(i));
    end


    fprintf("Con valores I = %i, J = %i, Método Integración %s\n",I,J,metodos{typeError})
    disp("")


    disp("Error 1")
    format shortEng
    disp(error_1)
    p = polyfit(log(h),log(error_1),1);

    disp("Ajuste 1")
    format("default")
    p1 = p(1);
    disp(p1);

    disp("Error 2")
    format shortEng
    disp(error_2)
    p = polyfit(log(h),log(error_2),1);

    disp("Ajuste 2")
    format("default")
    p2 = p(1);
    disp(p2);
end

