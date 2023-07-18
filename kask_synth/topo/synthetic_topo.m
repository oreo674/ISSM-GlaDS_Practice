function varargout = synthetic_topo(varargin)
% Compute synthetic topography
% outline = synthetic_topo(x) returns the outline
% [surf, bed] = synthetic_topo(x) returns surf and bed elevations

% CONSTANTS
ref_para = 0.05;
para = 0.05;
% x_max = 30e3;
x_max = 100e3;
y_scale = 3;
% z_scale = 1200/611;
z_scale = 3;
eps = 1e-10;

% DEFINE FUNCTIONS
surface = @(x, y) z_scale.*(100.*(x.*6e3/x_max + 200).^(1/4) + 1/60.*x*6e3/x_max - 100*(200).^(1/4)) + 1;
f = @(x, pa) (surface(x_max, 0) - pa*6e3)./x_max.^2 .* x.^2 + pa.*x*6e3/x_max;
g = @(y) 0.5e-6.*abs(y./y_scale).^3;
h = @(x, pa) ( (-4.5.*x/x_max + 5).*(surface(x, 0) - f(x, pa)))./(surface(x, 0) - f(x, ref_para) + eps);
bedfun = @(x, y, pa) f(x, pa) + g(y).*h(x, pa);
ginv = @(x) y_scale*(x./0.5e-6).^(1/3);
outline = @(x) ginv( (surface(x, 0) - f(x, ref_para))./(h(x, ref_para) + eps));

if nargin==1
    x = varargin{1};
    out = outline(x);
    varargout{1} = out;
else
    x = varargin{1};
    y = varargin{2};
    surf = surface(x, y);
    bed = bedfun(x, y, para);
    varargout{1} = surf;
    varargout{2} = bed;
end