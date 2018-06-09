% Program 6.6 Animation program for bridge using IVP solver
% Inputs: inter = [a b] time interval, 
%   ic = [y(1,1) y(1,2) y(1,3) y(1,4)],
% ic = y y' theta theta'
%   h = stepsize, p = steps per point plotted
% Calls a one-step method such as trapstep.m
% Example usage: tacoma2([0 1000],[0 0 0.001 0],.04,3);
% dito_656@hotmail.com
%ladotonia@gmail.com
function y_of_t = tacoma_rk(inter,ic,h,p)
clf                                  % clear figure window
h=inter(1);b=inter(2);n=ceil((b-h)/(h*p));% plot n points
y(1,:)=ic;                           % enter initial conds in y
t(1)=h;len=6;
figure(1)
set(gca,'XLim',[-8 8],'YLim',[-8 8], ...
   'XTick',[-8 0 8],'YTick',[-8 0 8], ...
   'Drawmode','fast','Visible','on','NextPlot','add');
cla;                                 % clear screen
axis square                          % make aspect ratio 1 - 1
road=line('color','b','LineStyle','-','LineWidth',5,...
    'erase','xor','xdata',[],'ydata',[]);
lcable=line('color','r','LineStyle','-','LineWidth',1,...
    'erase','xor','xdata',[],'ydata',[]);
rcable=line('color','r','LineStyle','-','LineWidth',1,...
    'erase','xor','xdata',[],'ydata',[]);
y_of_t = [];
for k=1:n
  for i=1:p
    t(i+1) = t(i)+h;
    y(i+1,:) = rk4step(t(i),y(i,:),h);
  end
  time(k)=t(p+1);
  y(1,:) = y(p+1,:);t(1)=t(p+1);
  z1(k)=y(1,1);z3(k)=y(1,3);
  c=len*cos(y(1,3));s=len*sin(y(1,3));
  set(road,'xdata',[-c c],'ydata',[-s-y(1,1) s-y(1,1)])
  set(lcable,'xdata',[-c -c],'ydata',[-s-y(1,1) 8])
  set(rcable,'xdata',[c c],'ydata',[s-y(1,1) 8])
  %drawnow; pause(h)
  y_of_t = [y_of_t; y(1,1),y(1,3)];
end
figure(2)
plot(time, y_of_t(:,1));
xlabel('Time');
ylabel('Angular Placement');
figure(3)
plot(time, y_of_t(:,2));

function y = trapstep(t,x,h)
%one step of the Trapezoid Method
z1=ydot(t,x);
g=x+h*z1;
z2=ydot(t+h,g);
y=x+h*(z1+z2)/2;

function y=rk4step(t,w,h)
%one step of the Runge-Kutta order 4 method
s1=ydot(t,w);
s2=ydot(t+h/2,w+h*s1/2);
s3=ydot(t+h/2,w+h*s2/2);
s4=ydot(t+h,w+h*s3);
y=w+h*(s1+2*s2+2*s3+s4)/6;

function ydot=ydot(t,y)
len=6; a=0.2; W=80; omega=2*pi*38/60;
a1=exp(a*(y(1)-len*sin(y(3))));
a2=exp(a*(y(1)+len*sin(y(3))));
ydot(1) = y(2);
ydot(2) = -0.01*y(2)-0.4*(a1+a2-2)/a+0.2*W*sin(omega*t);
ydot(3) = y(4);
ydot(4) = -0.01*y(4)+1.2*cos(y(3))*(a1-a2)/(len*a);