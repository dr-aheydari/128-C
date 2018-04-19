%% User Input

user_control_tf = input('Enter 1 to run Euler`s method with your own input, or 0 to continue to the output for Project 1: ');
if user_control_tf == 1
    a = input('what is the left endpoint (a)? ');
    b = input('what is the right endpoint (b)? ');
    y_a = input('What is y(a)? ');
    h = input('what is the step size (h)? ');
    y_prime = input('what is y_prime?(type @(t,y) [then the function]: ');

    if isempty(a)
        a = -5;
    end
    if isempty(b)
        b = 5;
    end
    if isempty(y_a)
        y_a = 1/(5^4 + 1);
    end
    if isempty(y_prime)
        y_prime = @(t,y) -4*(t^3)*(y^2);
    end
    if isempty(h)
        h = 10^-3;
    end

    % Runs Euler's method and stores the pairs of points (t_i, w_i) in the
    % matrix E, which is then split into vectors T and W for ease of plotting.
    E = euler_method(y_prime,a,b,y_a,h);
    T = E(:,1); W = E(:,2);
    
    figure
    plot(T,W);

else
%% Part (a): Euler's Method

    a = -5;
    b = 5;
    y_a = 1/(5^4 + 1);
    y_prime = @(t,y) -4*(t^3)*(y^2);
        
    E1 = euler_method(y_prime, a, b, y_a, 10^-3);
    T1 = E1(:,1); W1 = E1(:,2);
    % plot(T1,W1);
    % disp(E1);

    E2 = euler_method(y_prime, a, b, y_a, 10^-4);
    T2 = E2(:,1); W2 = E2(:,2);
    % plot(T2,W2);
    % disp(E2);


    E3 = euler_method(y_prime, a, b, y_a, 10^-5);
    T3 = E3(:,1); W3 = E3(:,2);

    figure
    plot(T1,W1, T2,W2, T3,W3);
    title('Euler`s Method Approximations');
    legend('h = 10^{-3}','h = 10^{-4}','h = 10^{-5}')
    % disp(E3);

    %% Modified Euler's Method

    EM1 = modi_euler(y_prime, a, b, y_a, 10^-3);
    TM1 = EM1(:,1); WM1 = EM1(:,2);
    % plot(TM1,WM1);
    % disp(EM1);

    EM2 = modi_euler(y_prime, a, b, y_a, 10^-4);
    TM2 = EM2(:,1); WM2 = EM2(:,2);
    % plot(TM2,WM2);
    % disp(EM2);

    EM3 = modi_euler(y_prime, a, b, y_a, 10^-5);
    TM3 = EM3(:,1); WM3 = EM3(:,2);
    % plot(TM3,WM3);
    % disp(EM3);
    
    % Note that you'll have to zoom in A LOT on the point (0,1) in order to see any
    % differences between the approximations in the following plot
    figure
    plot(TM1,WM1, TM2,WM2, TM3,WM3);
    title('Modified Euler`s Method Approximations');
    legend('h = 10^{-3}','h = 10^{-4}','h = 10^{-5}')

    %% Part (d): Richardson Extrapolation (Using results from Euler's Method)

    % This extracts the Euler's method approximations of y(0) from the
    % vectors of approximations computed in Part (a)
    results_from_part_a = [E1((length(E1)+1)/2,2) E2((length(E2)+1)/2,2) E3((length(E3)+1)/2,2)];
    
    vals = input('Input the first column (i.e., the results from Euler`s method) as a vector, or leave blank to use the values from Part (a): ');
    if isempty(vals)
        vals = results_from_part_a;
    end
    
    n = length(vals);
    A = zeros(n);
    A(:,1) = vals;

    % This algorithm computes the Richardson Extrapolation table using a
    % multiplier of 10 (each successive row divides h by 10) and stores it
    % in the matrix A. The first column of A are the Euler's method
    % approximations of y(0) for each value of h.
    for j = 2:n
        for i = j:n
            mult = 10^(j-1);
            A(i,j) = (mult * A(i,j-1) - A(i-1,j-1)) / (mult - 1);
        end
    end

    final_val = A(n,n);

    table_tf = input('Enter 1 to see the full table, or 0 to see only the final value: ');
    if table_tf
        disp('Richardson Extrapolation Table');
        disp(A);
    else
        disp(final_val);
    end

    error_tf = input('Enter 1 to see the table of errors, or 0 to see only the final error: ');
    if error_tf
        disp('Richardson Extrapolation Error Table');
        disp(tril(1 - A));
    else
        disp(1 - final_val);
    end
end
%% Functions

% This function computes the set of points for Euler's Method as described
% in Section 5.2 of our textbook.
function E = euler_method(f,a,b,y_a,h)
    N = floor((b-a)/h); % ceil( (b-a)/h ); % Compute N, rounding up to nearest integer if necessary
    t = a;
    w = y_a;
    E = zeros(N,2);
    E(1,:) = [t w];
    
    for i=1:N
        w = w + h*f(t,w); % Compute w_i
        t = a + i*h; % Compute t_i
        E(i+1,:) = [t w];
    end
end

% This function computes the set of points for Modified Euler's Method 
% as described on page 286 of our textbook.
function E = modi_euler(f,a,b,y_a,h)
    N = floor((b-a)/h); % ceil( (b-a)/h ); % Compute N, rounding up to nearest integer if necessary
    t = a;
    w = y_a;
    E = zeros(N,2);
    E(1,:) = [t w];
    
    for i=1:N
        w = w + (h/2)*( f(t,w) + f(a+i*h, w + h*f(t,w)) ); % Compute w_i
        t = a + i*h; % Compute t_i
        E(i+1,:) = [t w];
    end
end
