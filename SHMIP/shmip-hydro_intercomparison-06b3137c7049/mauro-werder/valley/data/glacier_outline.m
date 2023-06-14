function [x,y] = glacier_outline(n)

xend = 6e3;
theta = linspace(0,pi,n/2);
xx = xend/2 * (1-cos(theta));
y = [valley_outline(xx), -valley_outline(xx(end-1:-1:2))];
x = [xx, xx(end-1:-1:2)];