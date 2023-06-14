function p = model_para_valley(topopara, num)

if ~exist('num', 'var') || isempty(num)
    num = 100;
end
x0 = 0;
x1 = 6e3;

p = model_para(x0, x1, num);

%% physical
nil = p.coords * 0;  % zeros for all nodes

ev = 1e-3;
p.sigma = nil + ev*p.width/(p.rho_w*p.g);  % Volume storage capacity per unit length per unit pressure [m^2/Pa]

%% valley
p.topopara = topopara;

%find the surface
p.surf = @(x) 100.*(x+200).^(1/4) + 1/60.*x - 2e10.^(1/4) + 1;
%...the bed height
p.B  = (p.surf(6e3) - p.topopara.*6e3)/6e3.^2 * p.coords.^2 + p.topopara.*p.coords;
%... and the thickness
p.H = p.surf(p.coords) - p.B;
% width
p.width = real(2*(((p.H)./((-4.5.*p.coords/6e3 + 5) .* (p.H)./(p.H+eps())+eps())/0.5e-6).^(1/3)));

%% Initial conditions
p.S_IC = 3 + nil(1:end-1);
