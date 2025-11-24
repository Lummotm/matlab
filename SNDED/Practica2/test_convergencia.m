function [] = test_convergencia(I)
J = I+100;
N = I:10:J;
time_1 = zeros(1,length(N));
time_2 = zeros(1,length(N));
error_1 = zeros(1,length(N));
error_2 = zeros(1,length(N));

for i = 1:length(N)
    [~,~,time_1(i),error_1(i)] = els_finitos(1,4,N(i));
    [~,~,time_2(i),error_2(i)] = els_finitos(2,4,N(i));
end

fprintf("Con valores I = %i , J = %i",I,J)

error_1
error_2
end
