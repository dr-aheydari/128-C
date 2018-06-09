% Ali Heydari
% MAT 128 C
% Final Project, Computer Project1


%% part one: W = 80 ; checking for stability of the bridge
% Modified the given codes for user input

format long;
a = input('To run with default values, please enter 1, otherwise enter 2\n');

if (a == 1)
    tacoma([0 1000],[1 0 0.001 0],25000,5)


elseif ( a == 2)
       inter = input('Please enter time interval (example: [0 500])\n');
       ic = input ('Please enter initial condition (example: [1 0 0.001 0])\n');
       disp('To increase efficency, the number of steps will be calculated during function call\n');
       h = input('Please enter step size h (example: 0.05)');
       p = input('Please enter number of steps per point plotted (example: 5)\n'); 
       tacoma(inter,ic,h,p)
    
else
    
    disp(' You did not enter 1 or 2, please retry')
    
end

function tacoma(inter,ic,n,p)
clf                                  % clear figure window
h=(inter(2)-inter(1))/n;
y(1,:)=ic;                           % enter initial conds in y
t(1)=inter(1);len=6;
set(gca,'XLim',[-8 8],'YLim',[-8 8], ...
   'XTick',[-8 0 8],'YTick',[-8 0 8]);
cla;                                 % clear screen
axis square                          % make aspect ratio 1 - 1
road=animatedline('color','b','LineStyle','-','LineWidth',1);
lcable=animatedline('color','r','LineStyle','-','LineWidth',1);
rcable=animatedline('color','r','LineStyle','-','LineWidth',1);
for k=1:n
  for i=1:p
    t(i+1) = t(i)+h;
    y(i+1,:) = trapstep(t(i),y(i,:),h);
  end
  y(1,:) = y(p+1,:);t(1)=t(p+1);
  z1(k)=y(1,1);z3(k)=y(1,3);
  c=len*cos(y(1,3));s=len*sin(y(1,3));
  clearpoints(road);addpoints(road,[-c c],[-s-y(1,1) s-y(1,1)])
  clearpoints(lcable);addpoints(lcable,[-c -c],[-s-y(1,1) 8])
  clearpoints(rcable);addpoints(rcable,[c c],[s-y(1,1) 8])
  drawnow; pause(h)
end
end

function y = trapstep(t,x,h)
%one step of the Trapezoid Method
z1 = ydot(t,x);
g = x + h * z1;
z2 = ydot(t+h,g);
y = x + h * (z1 + z2) / 2;
end

function ydots = ydot(t,y)
len = 6; 
a = 0.2;
W=80; 
omega = 2 * pi * 38/60;
a1 = exp(a * (y(1) - len * sin(y(3))));
a2 = exp(a * (y(1) + len * sin(y(3))));
ydots(1) = y(2);
ydots(2) = -0.01 * y(2) - 0.4 * (a1 + a2 - 2) / ...
    a + 0.2 * W * sin(omega * t);
ydots(3) = y(4);
ydots(4) = -0.01 * y(4) + 1.2 * cos(y(3)) * (a1 - a2) / (len * a);

end