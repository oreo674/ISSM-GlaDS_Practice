function p = model_para_valley_e(run)

%%%%%%%%%%%%%% changing shape parameter
options=[0.05
         0
         -0.1
         -0.5
         -0.7];

topopara = options(run);

p = model_para_valley(topopara);

%% physical
nil = p.coords * 0;  % zeros for all nodes

p.M = nil +(1.158e-6*p.width); % source term given for all runs in section E... this amount is double run A6

%% Time
p.tspan = [0:10:500]*p.day;
