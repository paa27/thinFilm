%this fucntion needs the x vector, solution vector, h spacing of whole
%grid, and b
function [f,f3] = num_der(x,y,h,b)
%expand the vector with ones and the b values at the end (used 100 ghost points on each side)
o = ones(1,100);
bs = b*ones(1,100);
y = [o y(1:end) bs ];
%calculate the derivative using the gradient function and the spacing h
dydx = gradient(y,h);
dydx2 = gradient(dydx,h);
dydx3 = gradient(dydx2,h);
%resize vector back to orginal size
dydx = dydx(101:end-100);
dydx3= dydx3(101:end-100);
%plot(x,dydx)
%the derivative vector is passed back through the function in f
f = dydx;
f3 = dydx3;
end