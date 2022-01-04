close all
clear
clc

%Import data
stream=importdata("lidstream_re1000.txt");
vortex=importdata("lidvortex_re1000.txt");
par=importdata("lidparameters_re1000.txt");

stream=flipud(stream);
vortex=flipud(vortex);

%Extracting parameters
Nx=par.data(1);
Ny=par.data(2);
Lx=par.data(3);
Ly=par.data(4);
Px=par.data(5);
Py=par.data(6);
dt=par.data(7);
T=par.data(8);
Re=par.data(9);

%Grid spacings
dx=Lx/(Nx-1);
dy=Ly/(Ny-1);
halfx=0.5/dx+1;
halfy=0.5/dy+1;

%If no mid point, round up
if mod(halfx,1)~=0 || mod(halfy,1)~=0
    halfx=ceil(halfx);
    halfy=ceil(halfy);
end

%Defining grid
x=0:dx:Lx;
y=0:dy:Ly;

%Extracing data
ystream=stream(halfy,:);
xstream=stream(:,halfx);
xstream=xstream';

%Top and bottom boundaries
u(1)=(xstream(2)-xstream(1))/dy;
u(Ny)=(xstream(Ny)-xstream(Ny-1))/dy;

v(1)=(ystream(2)-ystream(1))/dy;
v(Ny)=(-ystream(Ny)+ystream(Ny-1))/dy;

%Using central difference scheme to find u and v
for i=2:Ny-1 
    u(i)=(-xstream(i-1)+xstream(i+1))/(2*dy); 
end

for i=2:Nx-1 
    v(i)=(-ystream(i-1)+ystream(i+1))/(2*dx); 
end

%Horizontal velocity u at y=0.5 at Reynolds
%numbers of 100, 400, 1000 and 3200 (all on the same plot please) 
%using a 161x161 grid and a domain
%size of Lx = Ly = 1;

figure(1)
plot(y,u,'-o','MarkerSize',2.5,'LineWidth', 1.5)
hold on
xlabel('y')
ylabel('u')
title('Velocity u with respect to y at x=0.5')
grid on
hold off

%Vertical velocity u at x=0.5 at Reynolds numbers of 100, 
%400, 1000 and 3200 (all on the same plot please) 
%using a 161x161 grid and a domain
%size of Lx = Ly = 1;
figure(2)
plot(x,v,'-o','MarkerSize',2.5,'LineWidth', 1.5)
hold on
xlabel('x')
ylabel('v')
title('Velocity v with respect to x at y=0.5')
grid on
hold off

%Contour plot of streamfunction at Re=100
figure(3)
[Ms,cs]=contourf(x,y,stream);
cs.LineWidth = 1.5;
cs.ShowText = 'on';
cs.LabelSpacing = 1000;
hold on
xlabel('x')
ylabel('y')
title('Contour plot of stream function')
grid on
hold off

%Contour plot of vorticity at Re=100
figure(4)
[Mv,cv]=contourf(x,y,vortex);
cv.LineWidth = 1.5;
cv.ShowText = 'on';
cv.LabelSpacing = 1000;
hold on
xlabel('x')
ylabel('y')
title('Contour plot of vorticity')
grid on
hold off

%Contour plot of streamfunction at (Lx, Ly) = (1, 2) and 
%(Lx, Ly) = (2, 1) at Re = 100
%figure(5)




%close all