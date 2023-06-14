function p = model_para_sqrt_c(run)

p = model_para_sqrt_b(5);

% options here
p.options = [0.25  % relative amp
             0.5
             1
             2];

%% physical
p.amp = p.options(run); % amplitude variations
p.M_steady = p.M;
p.M = @(t) double(moulin_source(t, p));

%% Time
p.tspan = [[-100:10:0], [1/24:1/24:30-1/24]]*p.day;
