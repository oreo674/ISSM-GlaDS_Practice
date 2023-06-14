function p = model_para_sqrt(num);
% para = model_para_sqrt(num);
%
% Defaults for sqrt Topo

if ~exist('num', 'var') || isempty(    num)
    num = 200;
end
x0 = 0;
x1 = 100e3;

p = model_para(x0, x1, num);

%% physical
nil = p.coords * 0;  % zeros for all nodes

p.M = 0; % is source term positive or negative
p.sigma = nil; % Volume storage capacity per unit length per unit pressure [m^2/Pa]


%% sqrt geometry
p.B = nil;
p.H = 6*( sqrt(p.coords+5e3) - sqrt(5e3) ) + 1;
p.width = 20e3;

%% Initial conditions
p.S_IC = 0.1 + nil(1:end-1);
