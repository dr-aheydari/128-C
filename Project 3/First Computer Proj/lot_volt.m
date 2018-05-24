% Ali Heydari
% Math 128 C
% HW3: CP 1
% group members: Alina Gataullina, Derek Mitchler


% timesteps
tsteps = [0 .01];

% col 1: rabbits, col 2: foxes
init_vak = [500 500]; 

% good ol' ODE solver
% [t, X] = ode45(@(t,X) lot_vot(t,X,0.01), tsteps, init_vak);


%% plotting the population vs time for multiple values of alpha

% plotting the phase plane

%% uncomment this to see the phase plots
% figure;
% plot(X(:,1), X(:,2));
% xlabel('Rabbit');
% ylabel('Fox');
% grid on;
% 
% title('ORIGINAL MODEL: Phase plane of Fox bs Rabbit when r_0 = 300, f_0 = 150');
%


% plot the population vs time
alphas = [ 5  10]; % Vector storing values of alphas
% 
  figure;
% 
 for i=1:length(alphas)
% 
[t, X] = ode45(@(t,X) lot_vot(t,X,alphas(i)), tsteps, init_vak);
  subplot(length(alphas),1,i);
% 
 plot(t, X(:,1), 'DisplayName', 'r(t)');
 hold on
 plot(t, X(:,2), 'DisplayName', 'f(t)')
xlabel('Time');
ylabel('Population');
legend('show')
title(['Using alpha = ' num2str(alphas(i)) ' in Pop. vs Time']); %')
% 
% xlabel('Rabbit');
% 
% ylabel('Fox');
% % how to do it : title(['y = sin(' num2str(i) 'x)'])
%  title(['Using alpha = ' num2str(alphas(i))]); %')
% 
% grid on;
% 
 end


% title('ORIGINAL MODEL: When r_0 = 300, f_0 = 150')

% 
% % this is for the first part where I compared the four alphas and 
% % - the phase plots 

% alphas = [0.01 2 5 10]; % Vector storing values of alphas
% 
%  figure;
% 
% for i=1:length(alphas)
% 
% 
% subplot(4,1,i);
% 
% plot(X(:,1), X(:,2));
% 
% xlabel('Rabbit');
% 
% ylabel('Fox');
% % how to do it : title(['y = sin(' num2str(i) 'x)'])
%  title(['Using alpha = ' num2str(alphas(i))]); %')
% 
% grid on;
% 
% end




% The given Lotka-Volterra model

function [x_prime] = lot_vot(t, X, alpha)

% OG MODEL
% 

% % initialize
x_prime = zeros(2,1);

% X(1) -> rabbit
x_prime(1) = 2*X(1) - alpha.*X(1).*X(2);

% X(2) -> fox
x_prime(2) = -X(2) + alpha.*X(1).*X(2);


% % NEW MODEL
% 
% 
% R = 400;
% 
% dXdt(1) = 2 * (1 - (X(1) / R)) * X(1) - alpha.*X(1).*X(2);
% 
% dXdt(2) = -X(2) + alpha.*X(1).*X(2);
% 


end
