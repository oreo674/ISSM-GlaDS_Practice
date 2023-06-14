function [bed, ice_thickness] = valley(x, y, bed_para)
% [bed, ice_thickness] = valley(x, y, bed_para)
%
% A valley shaped topography with fixed surface in which the bed can be changed (with
% bed_para).  It is setup such that changing bed_para does not change the outline.
%
% Default bed_para=300/6e3 leads to a glacier which mimics Bench-glacier's bed.
% Decreasing it will lead to thicker ice and eventually an overdeepened bed.

% Parameters

% domain length
xend = 6e3;

% surf para
beta = 1/4;
s1 = 100;
s2 = 100/xend;
sx0 = -200;
s3 = 1;

% bed para
g1 = .5e-6;
alpha = 3;
f20 = 300/xend;
if ~exist('bed_para', 'var') || isempty(bed_para)
    f2 = 300/xend;
else
    f2 = bed_para;
end

% surface
surf = s1*(x-sx0).^beta + s2*x - s1*(-sx0).^beta + s3;
s_xend = s1*(xend-sx0).^beta + s2*xend - s1*(-sx0).^beta + s3;

% bed:
h0 = -4.5.*(x./xend) + 5;
f10 = (s_xend - f20*xend)/xend^2;
f1 = (s_xend - f2*xend)/xend^2;

f0 = f10*x.^2 + f20*x;
f  = f1*x.^2  + f2*x;

h = h0 .* (surf - f)./(surf - f0);
h(isnan(h)) = 0;
h(isinf(h)) = 0;

bed = g1*abs(y).^alpha .* h + f;

% thickness
ice_thickness = surf-bed;
ice_thickness(ice_thickness<0) = 0; % round off errors
