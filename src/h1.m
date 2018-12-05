

load('/Users/Pablo/Documents/Repos/thinFilm/src/h0.mat','saved_sol')

m=0;
global b;
b = .1;
global xend;
xend = 20;

x = linspace(0,xend,200);
t = linspace(0,10,10);
global h0;
global h0prime;
global h0tripleprime;

h0 = saved_sol;

[h0prime,h0tripleprime] = num_der(x,h0,xend/200,b);

sol = pdepe(m,@pdefun,@icfun,@bcfun,x,t);

function [c,f,s] = pdefun(x,t,u,DuDx)
global h0;
global h0prime;
global h0tripleprime;
D=1.5;
q =0.4 ;
global b;
U=(1-b^3)/(1-b);
c=[0;1];
f = [zeros(length(h0)) ones(length(h0)); -h0.^3 (h0.^3.*D +q^2.*h0.^3) ]*DuDx + ... 
[zeros(length(h0)) zeros(length(h0)); ...
zeros(length(h0)) (-3.*h0.^2.*(h0tripleprime.^3 - ...
D.*h0prime)-(3.*h0.^2.*h0prime)+U.*ones(length(h0)))]*u;

s=[-1 0; q^2.*h0.^3 (-D*q^2.*h0.^3-q^4.*h0.^3) ].*u;
end

function u0 = icfun(x)
global b;
global xend;
u0 = zeros(length(x));
end

function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
global b;
pl = [0;0];
ql = [1;1];
pr = [0;0];
qr =[1;1];
end