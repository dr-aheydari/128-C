% Project 2
% MAT 128C
% Ali Heydari (team members Alina G., Derek M.)
%
% Runge-Kutta Fehlberg method


%  --------- this can be a user input program but for ease we just run it

% a = input('What is a?\n');
% 
% b = input('What is b?\n');
% y_a = input('y(a) = ?\n');
% h = input('Take a guess for h?\n');
% rtol = input('What is the tolerance (epsilon) ?\n'); 
% f = input('what is y_prime?(type @(t,y) [then the function]');

inp = input('press 1 for a = 5, press 2 for a = 10 \n');

if (inp == 1)
% first when a = 5;
 a = 0 ;
 b = 5;%-1 * a;
 y_a = 100 / ((a - 2)^4 + 1);
 h = .1; 
 rtol = 10 ^ (-11); 
 f = @(t,y) (-1 / 25) * ((t - 2)^(3)) * y^(2) ;

[w,t] = rkf(f,a,b,y_a,h,rtol);

end 

if (inp == 2)
    % Now when a = 10;
 a = -10 ;
 b = -1 * a;
 y_a = 100 / ((a - 2)^4 + 1);
 h = .1; 
 rtol = 10 ^ (-5); 
 f = @(t,y) (-1 / 25) * ((t - 2)^(3)) * y^(2) ;

[w,t] = rkf(f,a,b,y_a,h,rtol);
end

function [y t_array] = rkf(f,a,b,y_a,h,rtol)

E = zeros(63,3);

% bottom line minimum tolerance
last_tol = 1e-13;
const = .86;
% will use this later
check_pnt = 0;

% got the info from user
i = 1;
t_array(1) = a;
E(1,:) = [a ,y_a,h] 
t = a;

% Initial condition
y(1,:) = y_a;
w = y_a;
% our stop signal
FLAG = 0;

    while FLAG == 0

    % Checking if we are out of bounds or not
    % check if things are gucci
            if t + h > b
                   
                h = b - t;
                FLAG = 1;
            
            end
     
        k1 = h * f(t,w);
        k2 = h * f(t + ((1/4) * h), w + ((1/4) * k1));
        k3 = h * f(t + ((3/8) * h), w + ((3/32) * k1) + ((9/32)*k2));
        k4 = h * f(t + ((12/13) * h), w + ((1932/2197) * k1) - ...
            ((7200 / 2197) * k2) - ((7296 / 2197) * k3));
        k5 = h * f(t + h, w + ((439 / 216) * k1) - (8 * k2) + ...
            ((3680 / 513) * k3) - ((845 / 4104) * k4));
        k6 = h * f(t + ((1/2) * h), w - ((8 / 27) * k1) + (2 * k2) - ...
            ((3544 / 2565) * k3) + ((1859 / 4104) * k4) - ...
            ((11 / 40) * k5));
%    w = wi + h * (cw1 * s1 + cw3 * s3 + cw4 * s4 + cw5 * s5);

        w = w + (25 / 216) * k1 + (1408 / 2565) * k3 + (2197 / 4104) ...
            * k4 - (1 / 5) * k5;    
    
        z = w + h * ((16 / 135) * k1 + ((6656 /12825)  * k3) + ...
            ((28561/56430) * k4)+ (-9/50 * k5) + (2/55) * k6);


        R = (h) * norm((1/360) * k1 + (-128/4275) * k3 + (-2197 / 75240) ...
            * k4 + 1/50 * k5 + 2/55 * k6);
    
    % delta tolerance in the book
    delta = rtol * norm(w) + last_tol;
    % magic number that we vary h by
    q = const * (delta/R) ^ (1/5);
    
    if (R <= delta) % In case the tolerance is met
        t = t + h; 
        h = q * h;
        i = i + 1;
        t_array(i) = t;
        w = z;
        y(i,:) = z;
        check_pnt = 0;
            if (t > 0)
                E(i,:) = [t w h];
            end
    elseif (check_pnt == 0) % Tolerance is not met for the first time in this step
        
        h = const * h * (delta/R) ^ (1/7);
%         h = .1 * h;
        check_pnt = check_pnt + 1;
        FLAG = 0;
        
    else 
        
        h = 0.1 * h;
        FLAG = 0;
        
    end;
        
end;

fprintf("           t                   w                    h\n");
disp(E); 
end
