% Matlab script to illustrate change of coordinates from square to
% piecewise smooth annular domain for the noisy heteroclinic oscillator
% system.
%
% Note the code generates a trajectory moving clockwise.  For consistency
% with text the image is inverted to give CCW motion in figure 2.

close all
clear all

D=.01;
alpha=.01;
tmax=40;
dt=.01;
t=0:dt:tmax;
nt=length(t);
ntrace=5;
x1=nan(ntrace,nt); % cartesian coordinates x1,x2
x2=nan(ntrace,nt);
a1=nan(size(x1)); % transformed polaresque coordinates a1 (angle), a2 log radius
dw1=randn(size(x1));
dw2=randn(size(x2));
x1(:,1)=pi/4;
x2(:,1)=0;
a1(:,1)=0; % calculate a2 from trace later
f1=@(x1,x2)cos(x1).*sin(x2)+alpha*sin(2*x1);
f2=@(x1,x2)-sin(x1).*cos(x2)+alpha*sin(2*x2);

% some parameters
sdt=sqrt(2*D*dt);
pi2=pi/2;
epsilon=.01;

for j=1:nt
    xx1=x1(:,j)+dt*f1(x1(:,j),x2(:,j))+sdt*dw1(:,j);
    xx2=x2(:,j)+dt*f2(x1(:,j),x2(:,j))+sdt*dw2(:,j);
    xx1=(xx1<-pi2).*(-pi-xx1)+(xx1>pi2).*(pi-xx1)+(abs(xx1)<=pi2).*xx1;
    xx2=(xx2<-pi2).*(-pi-xx2)+(xx2>pi2).*(pi-xx2)+(abs(xx2)<=pi2).*xx2;
    x1(:,j+1)=xx1;
    x2(:,j+1)=xx2;
    da1=atan2(xx2,xx1)-atan2(x2(:,j),x1(:,j));
    da1=da1+2*pi*((da1<-pi)-(da1>pi));
    a1(:,j+1)=a1(:,j)+da1;
end
th=atan2(x2,x1);
absth=abs(th);
a2=(absth<=pi/4).*(x1-epsilon)./(pi2-epsilon)+...
    (th>pi/4).*(th<=3*pi/4).*(x2-epsilon)./(pi2-epsilon)+...
    (th<-pi/4).*(th>=-3*pi/4).*-(x2-epsilon)./(pi2-epsilon)+...
    (absth>3*pi/4).*-(x1-epsilon)./(pi2-epsilon);

a1cutoff=max(min(a1')); % all trajectories go this far.
a1shift=-2*pi*floor(a1cutoff/(2*pi)); % number of times around

%% figure 1b
figure(1)
clf
plot(x1',x2')
hold on
plot(x1(1,1),x2(1,1),'k.','MarkerSize',40)
set(line([-pi2 -pi2 pi2 pi2 -pi2],[pi2,-pi2,-pi2,pi2,pi2]),'Color','k','LineWidth',3)
set(line(epsilon*[-pi2 -pi2 pi2 pi2 -pi2],epsilon*[pi2,-pi2,-pi2,pi2,pi2]),'Color','k','LineWidth',1)
axis equal
axis off
shg
%print -dpdf het_osc_transform_figA

%% figure 1c
figure(2)
clf
subplot(4,1,1)
plot(a1'+a1shift,a2')
hold on
set(line([-4*pi a1shift+2*pi],[0 0]),'Color','k','LineWidth',3)
set(line([-4*pi a1shift+2*pi],[1 1]),'Color','k','LineWidth',3)
set(line([0 0],[0 1]),'Color','k','Linewidth',2,'LineStyle','--')
plot(a1(1,1)+a1shift,a2(1,1),'k.','MarkerSize',30)
axis off

