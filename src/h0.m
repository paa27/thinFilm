%ho solver
clear all
m=0;
global b;
b = .1;
global xend;
xend = 30;

x = linspace(0,xend,1000);
t = linspace(0,10,10);

options = odeset('AbsTol',1e-4);
sol = pdepe(m,@pdefun,@icfun,@bcfun,x,t,options);
saved_sol = sol(10,:,2);

save('/Users/Pablo/Documents/Repos/thinFilm/src/h0.mat','saved_sol')

for i=1:length(t)
    hold on
    plot(x,sol(i,:,2))
end 

function [c,f,s] = pdefun(x,t,u,DuDx)
D=1.5;
c=[0;1];
f=[0 1;-u(2).^3 u(2).^3*D]*DuDx+[0;-u(2).^3];
s=[-1;0].*u;
end

function u0 = icfun(x)
global b;
global xend;
u0 = [-(1-b)/(2*tanh(xend/2))*2.*(sech(-(x-xend/2))).^2*tanh(-(x-xend/2));
    (1-b)/(2*tanh(xend/2))*tanh(-(x-xend/2))+(1+b)/2];
end

function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
global b;
pl = [0;1];
ql = [1;1];
pr = [0;b.^3];
qr =[1;1];
end
