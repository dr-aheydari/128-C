% Ali Heydari
% MAT 128 C
% Final Project, Computer Project1


%% part Three: Finding W which gives magnifying factor of 100

% Modified the given codes

% n = 25000

format long

a = input('To run with default values, please enter 1, otherwise enter 2\n');

if (a == 1)
    % n = 25000 so h = 0.02
    tacoma_prt3([0 2500],[1 0 0.001 0],0.02,5);


elseif ( a == 2)
       inter = input('Please enter time interval (example: [0 500])\n');
       ic = input ('Please enter initial condition (example: [1 0 0.001 0])\n');
       n = input('Please enter the number of step (example: 2500)\n');
       p = input('Please enter number of steps per point plotted (example: 5)\n'); 
       h = (inter(2) - inter(1)) / (n * p);
       y_t = tacoma_prt3(inter,ic,h,p);
    
else
    
    disp(' You did not enter 1 or 2, please retry')
    
end


function y_t = tacoma_prt3(inter,ic,h,p)

a = inter(1);
b = inter(2);
n = ceil((b-a)/(h*p)); % getting the number of steps
count = 1;
MF = ones(1,10);
wind = ones(1,10);
for W = 50 : 100
    
    y_t = [];
    y = [];
    y(1,:) = ic;
    t = [];
    
    t(1) = a;
    for k = 1:n
        for i = 1:p
        t(i+1) = t(i)+h;
        y(i+1,:) = runge_kutta(t(i),y(i,:),h,W);
        end
    y(1,:) = y(p+1,:);
    t(1) = t(p+1);
    y_t = [y_t; y(1,1),y(1,3)];
    end
    mf = max(abs(y_t(:,2))) / abs(ic(1,3));
    fprintf(' for w = %i \n', W);
    fprintf(' Magnifying factor = %i\n', mf);
    
   wind(count) = W;
   MF(count) = mf;
   count = count+1;
   
end

hold on
    plot(wind,MF,'-o');
    title('W vs Magnifying Factor plot IC = [1,0,0.001,0]');
    xlabel('Wind Speed');
    ylabel('Magnifying Factor');
hold off

end

function y = runge_kutta(t,w,h,W)
%one step of the Runge-Kutta order 4 method
k1 = ydot(t,w,W);
k2 = ydot(t + h / 2, w + h * k1/2, W);
k3 = ydot(t + h / 2, w + h * k2/2, W);
k4 = ydot(t + h,w + h * k3, W);
y = w + h * (k1 + 2 * k2 + 2 * k3 + k4)/6;
end

function ydots=ydot(t,y,W)
len = 6;
a = 0.2;
omega = 2 * pi * 38 / 60;
a1 = exp(a * (y(1) - len * sin(y(3))));
a2 = exp(a * (y(1) + len * sin(y(3))));
ydots(1) = y(2);
ydots(2) = -0.02 * y(2) - 0.4 * (a1 + a2 - 2) / a + 0.2 ...
    * W * sin(omega * t);
ydots(3) = y(4);
ydots(4) = -0.02 * y(4) + 1.2 * cos(y(3)) * (a1 - a2) / (len * a);

end



