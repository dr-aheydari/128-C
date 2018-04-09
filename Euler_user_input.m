

a = input('what is a?\n');
b = input('what is b?\n');
y_a = input('y(a) = ?\n');
h = input('what is h?');
%f = @(t,y) exp(t-y);
f = input('what is y_prime?(type @(t,y) [then the function]');
euler_method(f,a,b,y_a,h)




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


