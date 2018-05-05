
% Setting all the values from the problem
 h = 0.001; % 10^-3
 M = 562000; % calculated for a = 5
 L = 1000; % a = 5
 a = 5;
 b = -5;
 
 N = (a - b) / h;

 err_bd = 0;

t;
counter = 0;
modulo = 100;

for i=0:N

    % finding t_i and error bound
    t = a + (i * h);
    err_bd = ((h*M)/(2*L)) * (exp(L*(t - a)) - 1);

    % to show every 100th value
    if mod(i, modulo) == 0
        fprintf(" error is %d\n", err_bd);
        counter = counter+1;
    end
end

fprintf( "the program showed only %d values\n", counter);
fprintf(" for every %d t_i 's\n", modulo);
