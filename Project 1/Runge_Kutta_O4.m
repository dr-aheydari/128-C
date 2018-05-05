% Runge-Kutta Order 4
% Alizei



a = input('what is a?\n');
b = input('what is b?\n');
y_a = input('y(a) = ?\n');
h = input('what is h?');
%f = @(t,y) exp(t-y);
f = input('what is y_prime?(type @(t,y) [then the function]');
modi_euler(f,a,b,y_a,h)






function E = modi_euler(f,a,b,y_a,h)
    N = ceil( (b-a)/h ); 
    % Compute N, rounding up to nearest integer if necessary
    
    t = a;
    w = y_a;
    E = zeros(N,2);
    E(1,:) = [t w];
    
    for i=1:N
        % here is the only change 
        % for now
        t = a + i*h; % Compute t_i
        %compute t_i+1
        q = a + (i-1)*h 
        
        K1=h*f(q,w);  
        K2 =h*f(q + h/2, w + K1 / 2);
        K3 =h*f(q + h/2, w + K2 /2); 
        K4 =h*f(q + h , w + K3);
        
        w = w + (1/6) * (K1 + (2 * K2)+(2 * K3)+K4); % compute w_i
        
        E(i+1,:) = [t w];
    end
end



