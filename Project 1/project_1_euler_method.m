%% Code for Example 1 in Section 5.2
% y_prime = @(t,y) y - t^2 +1;
% E = euler_method(y_prime, 0, 2, 0.5, 0.2);
% T = E(:,1); W = E(:,2);
% plot(T,W);
% disp(E);


%% Code for Computer Project 1
y_prime = @(t,y) 4*(t^3)*(y^2);

E1 = euler_method(y_prime, -5, 5, 1/(5^4 + 1), 10^-3);
T1 = E1(:,1); W1 = E1(:,2);
plot(T1,W1);
disp(E1);

% E2 = euler_method(y_prime, -5, 5, 1/(5^4 + 1), 10^-4);
% T2 = E2(:,1); W2 = E2(:,2);
% plot(T2,W2);
% disp(E2);


% E3 = euler_method(y_prime, -5, 5, 1/(5^4 + 1), 10^-5);
% T3 = E3(:,1); W3 = E3(:,2);
% plot(T3,W3);
% disp(E3);


% This function computes the set of points for Euler's Method as described
% in Section 5.2 of our textbook.
function E = euler_method(f,a,b,y_a,h)
    N = ceil( (b-a)/h ); % Compute N, rounding up to nearest integer if necessary
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
