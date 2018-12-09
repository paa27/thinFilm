load('/Users/Pablo/Documents/Repos/thinFilm/src/h0.mat','saved_sol')

m=0;
b = .1;
xend = 30;
q=0.2;
D=1.5;
tend=10;

x = linspace(0,xend,1000);
dx=xend/1000;
ts = linspace(0,tend,1e6);
dt = tend/1e6;

h0 = saved_sol;

[h0prime,h0tripleprime] = num_der(x,h0,dx,b);
h0 = [1 1 h0 b b];
h0prime(1)=0;
h0prime(end)=0;
h0prime = [0 0 h0prime 0 0];

g = zeros(length(ts),length(h0));
g(1,:) = h0prime;

for t=2:length(ts)
    der = zeros(1,length(h0));
    for x=4:length(h0)-3
        der(1,x) = h0(x).^3.*(g(t-1,x+2)-4.*g(t-1,x+1)+6.*g(t-1,x)-4.*g(t-1,x-1)+...
            g(t-1,x-2))*(dx).^(-4) + ...
            h0prime(x).^3.*(g(t-1,x+2)-2.*g(t-1,x+1)+2.*g(t-1,x-1)-...
        g(t-1,x-2)).*0.5.*(dx).^(-3)-...
        (D+2.*q.^2)*h0(x).^3.*(g(t-1,x+1)-2.*g(t-1,x)+g(t-1,x-1)).*(dx).^(-2)+...
        (2-h0prime(x).^3.*(D+q.^2)).*(g(t-1,x+1)-g(t-1,x-1))*(2.*dx).^(-1)+...
        h0(x).^3*q.^2.*(D+q.^2)*g(t-1,x);  
    end 
    g(t,:) = g(t-1,:)-dt.*der;
end
    
   
