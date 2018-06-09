

%% part b: Solving Burger's equation for k = 1/16, h = 0.01 and D = 0.2

% Modified the given codes for user input and option of showing graph or na

format long;
% initializing needed values

alf = 4;
bet = 3;
D = 0.2;
N_new = 16;
Graph = 1;
t = 1;
x = 1/2;
a = input('To run with default values, please enter 1, otherwise enter 2\n');

if (a == 1)
     
     w = burgers_1(Graph,0,1,0,1,100,16);
     
     
     figure(2);

p = 4:8;
% since we have te - tb = 1 - 0 = 1
N_new = arrayfun(@(x) 2^x, p);
k_new = arrayfun(@(x) x^-1,N_new);
% alf = 4, bet = 3 x = 1/2 t = 1;
err = zeros(5,1);
Graph  = 0;
w = burgers_1(Graph,0,1,0,1,100,k_new(1));

u = 2 * D * bet * pi * exp(-D * pi^2 * t) * sin(pi * x) / ...
    alf + bet*exp(-D*pi^2*t)*cos(pi*x);

for i = 1:5
    w = burgers_1(Graph,0,1,0,1,100,N_new(i));
    err(i) = abs(w(51, end) - u);
end

loglog(k_new, err, '-x', k_new, err)
xlabel('Log of Values of k')
ylabel('L/og Produced Error')
title('Log-log Error Plot for k = 2^{-p} for p = 4,5,6,7,8')

elseif ( a == 2)
       xl = input('Please enter left space interval (a.k.a xl)\n');
       xr = input ('Please enter right space interval (a.k.a xr)\n');
       tb = input('Please enter begining time(a.k.a tb): \n');
       te = input('Please enter end time (a.k.a te):\n'); 
       M = input('Please enter M:\n');
       N_new = input('Please enter N:\n');
       burgers_1(xl,xr,tb,te,M,N_new);
       
%        figure(2);

p = 4:8;
% since we have te - tb = 1 - 0 = 1
N_new = arrayfun(@(x) 2^x, p);
k_new = arrayfun(@(x) x^-1,k_new);
% alf = 4, bet = 3 x = 1/2 t = 1;
err = zeros(5,1);
Graph  = 0;
w = burgers_1(Graph,0,1,0,1,100,N_new(1));

u = 2 * D * bet * pi * exp(-D * pi^2 * t) * sin(pi * x) / ...
    alf + bet*exp(-D*pi^2*t)*cos(pi*x);

for i = 1:5
    u = act_fun(4,3,0.2,0.5,1);
    w = burgers_1(Graph,0,1,0,1,100,N_new(i));
    err(i) = abs(w(51, end) - u);
end

loglog(k_new, err, '-x', k_new, err)
xlabel('Log of Values of k')
ylabel('Log Produced Error')
title('Log-log Error Plot for k = 2^{-p} for p = 4,5,6,7,8')
    
else
    
    disp(' You did not enter 1 or 2, please retry')
    
end


function w = burgers_1(Graph,xl,xr,tb,te,M,N) 
alf = 4;
bet = 3;
D = .2; 
% f = @(x) 2 * D * bet * pi * sin(pi * x) ./ (alf + bet * cos(pi * x)); 
f=@(x) 2 * D * bet * pi * sin(pi * x) ./ (alf + bet * cos(pi * x)); 
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
if Graph  == 1
mesh(x,t,w');
xlabel('x');
ylabel('t');
title('Burgers equation solution for k = 1/16, h = 0.01, and D = 0.2');
else
    % do nothing
end
end

