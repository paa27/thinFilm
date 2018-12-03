function f = num_der(x,y,h,b)
%expand the vector with ones and the b values at the end (used 100 ghost points on each side)
o = ones(1,100);
bs = b*ones(1,100);
y = [o y(1:end) bs ];
%calculate the derivative using the gradient function and the spacing h
dydx = gradient(y,h);
%resize vector back to orginal size
dydx = dydx(101:end-100);
%plot(x,dydx)
%the derivative vector is passed back through the function in f
f = dydx;
end