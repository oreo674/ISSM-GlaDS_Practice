function y = valley_outline(x)
%  y = valley_outline(x)
%
% Gives the upper half of the glacier margin.

% Needs to be kept in sync with valley.m!

% domain length
xend = 6e3;

% surf para
beta = 1/4;
s1 = 100;
s2 = 100/xend;
sx0 = -200;
s3 = 0;

% bed para
g1 = .5e-6;
alpha = 3;
f20 = 300/xend;

surf = s1*(x-sx0).^beta + s2*x - s1*(-sx0).^beta + s3;
s_xend = s1*(xend-sx0).^beta + s2*xend - s1*(-sx0).^beta + s3;

% bed:
h0 = -4.5.*(x./xend) + 5;
f10 = (s_xend - f20*xend)/xend^2;

f0 = f10*x.^2 + f20*x;

y = (( (surf-f0)./h0  )/g1).^(1/alpha);