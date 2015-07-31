%Simplest Walker
%Last Modified 3/18/2009
% Josh Petersen
% Based on rimlesswheel.m by Pranav Bhounsule

% http://ruina.tam.cornell.edu/research/topics/locomotion_and_robotics/ranger/ranger_paper/Reports/Ranger_Robot/control/simulator/simplest_walker.m
function simplest_walker

clc
clear all
close all
format long


M = 1000; m = 1.0; I = 0.00; l = 1.0; w = 0.0; 
c = 1.0;  r = 0.0; d = 0.00; g = 1.0; gam = 0.009; 

steps=50;

GL_DIM = [M m c I g l w r d gam];

q10=.2; u10=-.2;
q20=pi; u20=0;
q30=.4; u30=-.4;
q40=pi; u40=0;

z0=[q10 u10 q20 u20 q30 u30 q40 u40];


%%%%% Root finding, Period one gait %%%%
options = optimset('TolFun',1e-6,'TolX',1e-6,'Display','off');
[zstar,fval,exitflag,output,jacob] = fsolve(@fixedpt,z0,options,GL_DIM);
if exitflag == 1
    disp('Fixed points are');
    zstar
else
    error('Root finder not converged, change guess or change system parameters')
end

%%%% Stability, using linearised eigenvalue %%%
J=partialder(@onestep,zstar,GL_DIM);
disp('EigenValues for linearized map are');
eig(J)

[z,t]=onestep(zstar,GL_DIM,steps);

q1 = z(:,1);
q2 = z(:,3); 
q3 = z(:,5);
q4 = z(:,7);


%Plotting
figure(1);
axis([-2 2 -2 2]);

for i=1:length(t)

xh=l*sin(q1(i)+q2(i))-d*sin(q1(i))-r*q1(i);
yh=-l*cos(q1(i)+q2(i))+d*cos(q1(i))+r;
xa1=xh-sin(q1(i)+q2(i));
ya1=yh+cos(q1(i)+q2(i));

xa2=xh+l*sin(q3(i)-q1(i)-q2(i));
ya2=yh+l*cos(q3(i)-q1(i)-q2(i));
figure(1)
plot([xh xa1],[yh ya1]);%plots the right half of the walker
hold on;
plot([xh xa2],[yh ya2]);% plots the left half
plot([xh-1 xh+1],[0 0],'black')%plots the ground line
%plot(xh,yh,'ko-','LineWidth',3)

axis([xh-1 xh+1 -1 1.5]);
axis off;
hold off;

end

figure(2)
hold on
plot(t,q1,'r');
plot(t,q3,'b');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('Stance Angle','Swing Angle')
hold off

%%%%Functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zdiff=fixedpt(z0,GL_DIM)
zdiff=onestep(z0,GL_DIM)-z0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J=partialder(FUN,z,GL_DIM)
pert=1e-5;
 
%%% Using central difference, accuracy quadratic %%%
for i=1:length(z)
    ztemp1=z; ztemp2=z;
    ztemp1(i)=ztemp1(i)+pert; 
    ztemp2(i)=ztemp2(i)-pert; 
    J(:,i)=(feval(FUN,ztemp1,GL_DIM)-feval(FUN,ztemp2,GL_DIM)) ;
end
J=J/(2*pert);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z,t]=onestep(z0,GL_DIM,steps)

M = GL_DIM(1);  m = GL_DIM(2); c = GL_DIM(3);   
I = GL_DIM(4);  g = GL_DIM(5); l = GL_DIM(6);   
w = GL_DIM(7);  r = GL_DIM(8); d = GL_DIM(9);   
gam = GL_DIM(10);

flag = 1;
if nargin<2
    error('need more inputs to onestep');
elseif nargin<3
    flag = 0; %send only last state
    steps = 1;
end

t0 = 0; 
dt = 5;
t_ode = t0;
z_ode = z0;

for i=1:steps
    options=odeset('abstol',2.25*1e-14,'reltol',2.25*1e-14,'events',@collision);
    tspan = linspace(t0,t0+dt,1000);
    [t_temp, z_temp, tfinal] = ode113(@ranger_ss_simplest,tspan,z0,options,GL_DIM);

    zplus=heelstrike_ss_simplest(t_temp(end),z_temp(end,:),GL_DIM);
    
    z0 = zplus;
    t0 = t_temp(end);
    
    %%% dont include the first point
    t_ode = [t_ode; t_temp(2:end); t0];
    z_ode = [z_ode; z_temp(2:end,:); z0];
    
end

z = [zplus(1:2) pi 0 zplus(5:6) pi 0];

if flag==1
   z=z_ode;
   t=t_ode;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function zdot = ranger_ss_simplest(t,z,GL_DIM);

q1=z(1);
u1=z(2);
q2=z(3);
u2=z(4);
q3=z(5);
u3=z(6);
q4=z(7);
u4=z(8);

Ta=0;
Th=0;
Thip=0;

M = GL_DIM(1);  m = GL_DIM(2); c = GL_DIM(3);   
I = GL_DIM(4);  g = GL_DIM(5); l = GL_DIM(6);   
w = GL_DIM(7);  r = GL_DIM(8); d = GL_DIM(9);   
gam = GL_DIM(10);                                       

M11 = 2*M*d*l*cos(q2)-2*I-2*m*l^2-2*m*w^2-2*m*d*c*cos(-q3+q2)+2*m*d*w*sin(-q3+q2)+2*m*l*c*cos(q3)+2*m*l*w*sin(q3)-M*r^2-2*m*r^2-M*l^2-M*d^2-2*m*d^2+2*m*l*c-2*m*r*c*cos(-q3+q1+q2)+2*m*r*w*sin(-q3+q1+q2)-4*m*r*d*cos(q1)+4*m*r*l*cos(q1+q2)-2*m*r*c*cos(q1+q2)+2*m*r*w*sin(q1+q2)-2*M*r*d*cos(q1)+2*M*r*l*cos(q1+q2)+4*m*d*l*cos(q2)-2*m*d*c*cos(q2)+2*m*d*w*sin(q2)-2*m*c^2; 
M13 = I+m*r*c*cos(-q3+q1+q2)-m*d*w*sin(-q3+q2)-m*l*w*sin(q3)-m*l*c*cos(q3)-m*r*w*sin(-q3+q1+q2)+m*w^2+m*c^2+m*d*c*cos(-q3+q2); 

M31 = -I+m*l*w*sin(q3)+m*l*c*cos(q3)-m*c^2-m*w^2-m*c*d*cos(-q3+q2)+m*w*d*sin(-q3+q2)-m*c*cos(-q3+q1+q2)*r+m*w*sin(-q3+q1+q2)*r; 
M33 = m*w^2+m*c^2+I; 

RHS1 = -m*g*((d*sin(q1)-l*sin(q1+q2)+c*sin(q1+q2)+w*cos(q1+q2))*cos(gam)-(r+d*cos(q1)-l*cos(q1+q2)+c*cos(q1+q2)-w*sin(q1+q2))*sin(gam))-m*g*((d*sin(q1)-l*sin(q1+q2)+c*sin(-q3+q1+q2)+w*cos(-q3+q1+q2))*cos(gam)-(r+d*cos(q1)-l*cos(q1+q2)+c*cos(-q3+q1+q2)-w*sin(-q3+q1+q2))*sin(gam))-M*g*((d*sin(q1)-l*sin(q1+q2))*cos(gam)-(r+d*cos(q1)-l*cos(q1+q2))*sin(gam))+m*((-d*sin(q1)+l*sin(q1+q2)-c*sin(q1+q2)-w*cos(q1+q2))*(((-d*cos(q1)+l*cos(q1+q2))*u1+l*cos(q1+q2)*u2)*u1+(l*cos(q1+q2)*u1+l*cos(q1+q2)*u2)*u2-(u1+u2)^2*(c*cos(q1+q2)-w*sin(q1+q2)))-(r+d*cos(q1)-l*cos(q1+q2)+c*cos(q1+q2)-w*sin(q1+q2))*(((d*sin(q1)-l*sin(q1+q2))*u1-l*sin(q1+q2)*u2)*u1+(-l*sin(q1+q2)*u1-l*sin(q1+q2)*u2)*u2-(u1+u2)^2*(-c*sin(q1+q2)-w*cos(q1+q2))))+m*((-d*sin(q1)+l*sin(q1+q2)-c*sin(-q3+q1+q2)-w*cos(-q3+q1+q2))*(((-d*cos(q1)+l*cos(q1+q2))*u1+l*cos(q1+q2)*u2)*u1+(l*cos(q1+q2)*u1+l*cos(q1+q2)*u2)*u2-(-u3+u1+u2)^2*(c*cos(-q3+q1+q2)-w*sin(-q3+q1+q2)))-(r+d*cos(q1)-l*cos(q1+q2)+c*cos(-q3+q1+q2)-w*sin(-q3+q1+q2))*(((d*sin(q1)-l*sin(q1+q2))*u1-l*sin(q1+q2)*u2)*u1+(-l*sin(q1+q2)*u1-l*sin(q1+q2)*u2)*u2-(-u3+u1+u2)^2*(-c*sin(-q3+q1+q2)-w*cos(-q3+q1+q2))))+M*((-d*sin(q1)+l*sin(q1+q2))*(((-d*cos(q1)+l*cos(q1+q2))*u1+l*cos(q1+q2)*u2)*u1+(l*cos(q1+q2)*u1+l*cos(q1+q2)*u2)*u2)-(r+d*cos(q1)-l*cos(q1+q2))*(((d*sin(q1)-l*sin(q1+q2))*u1-l*sin(q1+q2)*u2)*u1+(-l*sin(q1+q2)*u1-l*sin(q1+q2)*u2)*u2)); 
RHS3 = -m*g*c*sin(-gam-q3+q1+q2)-m*g*w*cos(-gam-q3+q1+q2)-Th+Thip-m*w*l*u2^2*cos(q3)-2*m*w*u1*l*u2*cos(q3)+m*c*l*u1^2*sin(q3)+m*c*u1^2*d*sin(-q3+q2)-m*w*l*u1^2*cos(q3)+m*c*l*u2^2*sin(q3)+m*w*u1^2*d*cos(-q3+q2)+2*m*c*u1*l*u2*sin(q3); 

MM = [M11 M13; M31 M33];                             

RHS = [RHS1; RHS3];                      

X = MM \ RHS;                                    

ud1 = X(1);                                       
ud2 = 0;                                       
ud3 = X(2);                                       
ud4 = 0;                                        

zdot = [u1 ud1 u2 ud2 u3 ud3 u4 ud4]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function zplus=heelstrike_ss_simplest(t,z,GL_DIM)      

r1 = z(1);   v1 = z(2);                         
r2 = z(3);   v2 = z(4);                         
r3 = z(5);   v3 = z(6);                         
r4 = z(7);   v4 = z(8);                       

q1 = r1 + r2 - r3 - r4;                         
q2 = r4;                                        
q3 = -r3;                                       
q4 = r2; % DEPENDS                           

M = GL_DIM(1);  m = GL_DIM(2); c = GL_DIM(3);   
I = GL_DIM(4);  g = GL_DIM(5); l = GL_DIM(6);   
w = GL_DIM(7);  r = GL_DIM(8); d = GL_DIM(9);   
gam = GL_DIM(10);                                                      

M11 = 2*c^2*m-2*m*l*w*sin(q3)-2*m*d*w*sin(-q3+q2)-2*m*l*c*cos(q3)+M*l^2+2*I+2*w^2*m+2*l^2*m-4*m*d*l*cos(q2)+2*m*d*c*cos(q2)-2*m*d*w*sin(q2)+M*d^2+2*m*d^2+2*m*r^2-2*c*l*m-2*M*d*l*cos(q2)+M*r^2-4*m*r*l*cos(q1+q2)+4*m*r*d*cos(q1)+2*m*r*c*cos(q1+q2)-2*m*r*w*sin(q1+q2)+2*m*r*c*cos(-q3+q1+q2)-2*m*r*w*sin(-q3+q1+q2)-2*M*r*l*cos(q1+q2)+2*M*r*d*cos(q1)+2*m*d*c*cos(-q3+q2); 
M13 = m*d*w*sin(-q3+q2)+m*l*c*cos(q3)+m*l*w*sin(q3)+m*r*w*sin(-q3+q1+q2)-m*r*c*cos(-q3+q1+q2)-m*d*c*cos(-q3+q2)-c^2*m-w^2*m-I; 

M31 = -m*l*w*sin(q3)+m*c^2+m*w^2-m*l*c*cos(q3)+m*c*d*cos(-q3+q2)-m*w*d*sin(-q3+q2)+m*c*cos(-q3+q1+q2)*r-m*w*sin(-q3+q1+q2)*r+I; 
M33 = -m*c^2-m*w^2-I; 

RHS1 = m*((-d*sin(-r3+r1+r2-r4)+l*sin(-r3+r1+r2)-c*sin(-r3+r1+r2)-w*cos(-r3+r1+r2))*((-d*sin(r1)+l*sin(r1+r2))*v1+l*sin(r1+r2)*v2+(-v3+v1+v2)*(-c*sin(-r3+r1+r2)-w*cos(-r3+r1+r2)))-(r+d*cos(-r3+r1+r2-r4)-l*cos(-r3+r1+r2)+c*cos(-r3+r1+r2)-w*sin(-r3+r1+r2))*((l*cos(r1+r2)-d*cos(r1)-r)*v1+l*cos(r1+r2)*v2-(-v3+v1+v2)*(c*cos(-r3+r1+r2)-w*sin(-r3+r1+r2))))+m*((-d*sin(-r3+r1+r2-r4)+l*sin(-r3+r1+r2)-c*sin(r1+r2)-w*cos(r1+r2))*((-d*sin(r1)+l*sin(r1+r2))*v1+l*sin(r1+r2)*v2+(v1+v2)*(-c*sin(r1+r2)-w*cos(r1+r2)))-(r+d*cos(-r3+r1+r2-r4)-l*cos(-r3+r1+r2)+c*cos(r1+r2)-w*sin(r1+r2))*((l*cos(r1+r2)-d*cos(r1)-r)*v1+l*cos(r1+r2)*v2-(v1+v2)*(c*cos(r1+r2)-w*sin(r1+r2))))+M*((-d*sin(-r3+r1+r2-r4)+l*sin(-r3+r1+r2))*((-d*sin(r1)+l*sin(r1+r2))*v1+l*sin(r1+r2)*v2)-(r+d*cos(-r3+r1+r2-r4)-l*cos(-r3+r1+r2))*((l*cos(r1+r2)-d*cos(r1)-r)*v1+l*cos(r1+r2)*v2))+I*(2*v1+2*v2-v3); 
RHS3 = m*c*v1*d*cos(r2)-m*w*v1*d*sin(r2)-m*l*v1*c-m*l*v2*c+m*c^2*v2+m*c^2*v1+m*w^2*v1+m*w^2*v2+m*r*v1*c*cos(r1+r2)-m*r*v1*w*sin(r1+r2)+I*v1+I*v2; 

MM = [M11 M13; M31 M33];                             

RHS = [RHS1; RHS3];                      

X = MM \ RHS;                                    

u1 = X(1);                                       
u2 = 0;                                       
u3 = X(2);                                       
u4 = 0;       %DEPENDS                       

zplus = [q1 u1 q2 u2 q3 u3 q4 u4]; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gstop, isterminal,direction]=collision(t,z,GL_DIM)

M = GL_DIM(1);  m = GL_DIM(2); c = GL_DIM(3);   
I = GL_DIM(4);  g = GL_DIM(5); l = GL_DIM(6);   
w = GL_DIM(7);  r = GL_DIM(8); d = GL_DIM(9);   
gam = GL_DIM(10);   

q1 = z(1); q2 = z(3); 
q3 = z(5); q4 = z(7); 

gstop = -l*cos(q1+q2)+d*cos(q1)+l*cos(-q3+q1+q2)-d*cos(q1+q2-q3-q4);
if (q3>-0.05) %no collision detection for foot scuffing
    isterminal = 0;
else
    isterminal=1; %Ode should terminate is conveyed by 1, if you put 0 it goes till the final time u specify
end
direction=-1; % The t_final can be approached by any direction is indicated by this