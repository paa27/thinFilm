%ho solver
m=0;
x = linspace(0,10,10000);
t = [0 0.005 0.01 0.05 0.1 0.5 1 1.5 2];

sol = pdepe(m,@pdefun,@icfun,@bcfun,x,t);

for i=1:length(t)
    hold on
    plot(x,sol(i,:,2))
end 

function [c,f,s] = pdefun(x,t,u,DuDx)
D=5;
c=[0;1];
f=[-1 1;-u(2).^3 u(2).^3*D]*DuDx+[0;-u(2).^3];
s=[0;0];
end

function u0 = icfun(x)
b=0.1;
xend=10;
u0 = [-(1-b)/(2*tanh(xend/2))*2.*(sech(-(x-xend/2))).^2*tanh(-(x-xend/2));
    (1-b)/(2*tanh(xend/2))*tanh(-(x-xend/2))+(1+b)/2];
end

function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
b=0.1;
pl = [0;1];
ql = [1;1];
pr = [0;b.^3];
qr =[1;1];
end

