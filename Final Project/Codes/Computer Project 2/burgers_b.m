

%% part b: Solving Burger's equation for k = h = 0.02

% Modified the given codes for user input

format long;

a = input('To run with default values, please enter 1, otherwise enter 2\n');

if (a == 1)
   
     burgers_1(0,1,0,1,50,50);

elseif ( a == 2)
       xl = input('Please enter left space interval (a.k.a xl)\n');
       xr = input ('Please enter right space interval (a.k.a xr)\n');
       tb = input('Please enter begining time(a.k.a tb): \n');
       te = input('Please enter end time (a.k.a te):\n'); 
       M = input('Please enter M:\n');
       N = input('Please enter N:\n');
       burgers_1(xl,xr,tb,te,M,N);
    
else
    
    disp(' You did not enter 1 or 2, please retry')
    
end

% Program 8.7 Implicit Newton solver for Burgers equation 
% input: space interval [xl,xr], time interval [tb,te], 
%        number of space steps M, number of time steps N 
% output: solution w 
% Example usage: w=burgers(0,1,0,2,20,40) 

function w = burgers_1(xl,xr,tb,te,M,N) 
alf = 5;
bet = 4;
D = .2; 
% f = @(x) 2 * D * bet * pi * sin(pi * x) ./ (alf + bet * cos(pi * x)); 
f = @(x)sin(2* pi * x);

l = @(t) 0 * t; 
r = @(t) 0 * t; 
 h = (xr - xl) / M;
 k = (te - tb) / N;
m = M + 1; 
n = N; 
sigma = D * k / (h * h);
w(:,1) = f(xl + (0:M) * h)'; % initial conditions
w1=w; 

for j=1:n 
  for it = 1:3 % Newton iteration 
    DF1 = zeros(m,m);
    DF2 = zeros(m,m); 
    
    DF1 = diag(1 + 2 * sigma * ones(m,1)) + diag(-sigma * ones(m - 1,1),1); 
    DF1 = DF1 + diag(-sigma * ones(m - 1,1),-1); 
    DF2 = diag([0;k * w1(3:m) / (2 * h);0]) - diag([0;k * ...
        w1(1:(m-2)) / (2 * h);0]); 
    DF2 = DF2 + diag([0;k * w1(2:m - 1) / (2 * h)],1) - diag([k * ...
        w1(2:m - 1) / (2 * h);0],-1); 
    DF = DF1 + DF2; 
    F = -w(:,j) + (DF1 + DF2 / 2) * w1; % Using Lemma 8.11 
    DF(1,:) = [1 zeros(1,m - 1)]; % Dirichlet conditions for DF 
    DF(m,:) = [zeros(1,m - 1) 1]; 
    F(1) = w1(1) - l(j);
    F(m) = w1(m) - r(j); % Dirichlet conditions for F 
    w1 = w1 - DF \ F; 
  end 
  w(:,j+1) = w1; 
end 

x = xl + (0:M) * h;
t = tb + (0:n) * k; 
mesh(x,t,w');
xlabel('x');
ylabel('t');
zlabel('Function Values');
title('Burgers equation solution for k = h = 0.02, and D = 0.2');

end


