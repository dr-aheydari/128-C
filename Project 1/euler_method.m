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
