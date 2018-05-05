% Project 2
% MAT 128C
% Ali Heydari (team members Alina G., Derek M.)
%
% Runge-Kutta Fehlberg method
% this could be an interactive code but for ease of use, we have made it a
% menu option

% --------------------------------------------------------------------- %
% format long;
% a = input('What is a?\n');
% b = input('What is b?\n');
% y_a = input('y(a) = ?\n');
% h = input('Take a guess for h?\n');
% tol = input('What is the tolerance (epsilon) ?\n'); 
% %f = @(t,y) exp(t-y);
% f = input('what is y_prime?(type @(t,y) [then the function])');
% --------------------------------------------------------------------- %

inp = input('Please press 1 to see the graph for a == 5 \n or press 2 to see for a == 10 \n');

% when a == 5
if (inp == 1)
% this just to initilize t
a = -5; % actual a value is a = 5
b = 5;
% given
y_a = 100 / ((a - 2)^4 + 1);
% initial guess
hmax = .25;
% our pick for smallest h step
hmin = .0001;
%first tolerance
tol = 10 ^ (-7); 
 f = @(t,y) (-1 / 25) * ((t - 2)^(3)) * y^(2) ;

[E] = RKF(f,a,b,y_a,hmax, hmin,tol);


plot(E(:,1),E(:,2));
title('RKF approximation');
xlabel('Time');
ylabel('Approximation using RK-Fehlberg');

hold on
plot(E(:,1),E(:,3));

tol = 10 ^ (-11); 

[E] = RKF(f,a,b,y_a,hmax, hmin,tol);

 plot(E(:,1),E(:,2));
% title('RKF approximation');
% xlabel('Time');
% ylabel('Approximation using RK-Fehlberg');
 plot(E(:,1),E(:,3));

 fplot(@(x) 100 / ((x - 2)^4 + 1),[-5 5],'r')
grid on

legend('esp = 10^-7','h steps for 10^-7', 'eps = 10^-11', ' h steps for 10^-11', 'Exact Solution');

hold off

end 

% when a  == 10

if (inp == 2)

a = -10; 
b = 10;
% given
y_a = 100 / ((a - 2)^4 + 1);
% initial guess
hmax = .1;
% our pick for smallest h step
hmin = .00001;
%first tolerance
tol = 10 ^ (-7); 
f = @(t,y) (-1 / 25) * ((t - 2)^(3)) * y^(2) ;

[E] = RKF(f,a,b,y_a,hmax, hmin,tol);


plot(E(:,1),E(:,2));
title('RKF approximation');
xlabel('Time');
ylabel('Approximation using RK-Fehlberg');
hold on
plot(E(:,1),E(:,3));

tol = 10 ^ (-11); 
[E] = RKF(f,a,b,y_a,hmax, hmin,tol);

plot(E(:,1),E(:,2));
title('RKF approximation');
xlabel('Time');
ylabel('Approximation using RK-Fehlberg');
plot(E(:,1),E(:,3));

fplot(@(x) 100 / ((x - 2)^4 + 1),[-10 10],'r')

grid on

legend('esp = 10^-7','h steps for 10^-7', 'eps = 10^-11', ' h steps for 10^-11', 'Exact Solution');

hold off
end 



function [E] = RKF(f,a,b,y_a, hmax, hmin,tol)
% we got the user input in the main

t = a;
w = y_a;

% something
h = hmax;
% for the while loop
FLAG = 1;
% array of points
E = zeros(10,3);
E(1,:) = [t w h];
 i = 1;
% if we need to do calculations
    while (FLAG == 1)
        
       
        % RKF Algortihm
        %disp(h);
        %disp("Break Break Break");
        k1 = h * f(t,w);
%         disp(k1);
%         disp("BREAK BREAK BREAK")
        k2 = h * f(t + ((1/4) * h), w + ((1/4) * k1));
%         disp(k2);
%         disp("BREAK BREAK BREAK")
        k3 = h * f(t + ((3/8) * h), w + ((3/32) * k1) + ((9/32)*k2));
%         disp(k3);
%         disp("BREAK BREAK BREAK")
        k4 = h * f(t + ((12/13) * h), w + ((1932 / 2197) * k1) - ...
            ((7200 / 2197) * k2) + ((7296 / 2197) * k3));
%         disp(k4);
%         disp("BREAK BREAK BREAK")
        k5 = h * f(t + h, w + ((439 / 216) * k1) - (8 * k2) + ...
            ((3680 / 513) * k3) - ((845 / 4104) * k4));
%         disp(k5);
%         disp("BREAK BREAK BREAK")
        k6 = h * f(t + ((1/2) * h), w - ((8 / 27) * k1) + (2 * k2) - ...
            ((3544 / 2565) * k3) + ((1859 / 4104) * k4) - ...
            ((11 / 40) * k5));
%         disp(k6);
%         disp("BREAK BREAK BREAK")
        
        R = (1/h) * abs(  ((1/360) * k1) - ((128 / 4275) * k3) - ...
            ((2197 / 75240) * k4) + ((1 / 50) * k5) + ((2 / 55) * k6));
        
%         disp(R);
        
        if (R <= tol)
%             disp(h);
            t = t + h; % approximation accepted
            w = w + (((25 / 216) * k1) + ((1408 / 2565) * k3) ...
                + ((2197 / 4104) * k4) - ((1 / 5) * k5));
        % OUTPUT HERE ACCORDING TO THE BOOK
         E(i+1,:) = [t w h];
%          FLAG = 1;
        end
        
         % in the book   
        q = 0.84 * (tol / R)^(1 / 4);
        
        if (q <= 0.1)
            
            h = 0.1 * h;
            
        elseif (q >= 4)
            
            h = 4 * h; 
            
        else
            
            h = q * h;
            
        end
    
         if (h > hmax)
             
             h = hmax;
         
         end
         
         
         
           if ( t >= b )
            
            FLAG = 0;
            
           elseif (t + h > b)
           
               h = b - t;
       
            elseif (h < hmin)
             
            FLAG = 0;
            fprintf(" min h exceeded ! \n");
            
            end
%         if (R > tol)
%            
%             fprintf("THIS IS WAY TOO BIG %i instead of h \n", h);
%             FLAG = 0;
%         end
%         
    
        E(i+1,:) = [t w h];
        
        i = i + 1;
        OUT = E;
    end
    
%     disp("Columns are : t, w, h respectively");
    end
    
