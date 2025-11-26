I = 10;
J = I+100;
N = I:10:J;
h = 1./N;

time_1 = zeros(1,length(N));
time_2 = zeros(1,length(N));
error_1 = zeros(1,length(N));
error_2 = zeros(1,length(N));

for j = 1:4

    typeError = j;
    for i = 1:length(N)
        [~,~,time_1(i),error_1(i)] = els_finitos(1,typeError,N(i));
        [~,~,time_2(i),error_2(i)] = els_finitos(2,typeError,N(i));
    end


    fprintf("Con valores I = %i , J = %i , TypeError = %i\n",I,J,typeError)

    % format shortEng
    % disp(error_1)
    % disp(error_2)
    % format

    p = polyfit(log(h),log(error_1),1);
    p1 = p(1);
    p = polyfit(log(h),log(error_2),1);
    p2 = p(1);

    disp(p1);
    disp(p2);
end

