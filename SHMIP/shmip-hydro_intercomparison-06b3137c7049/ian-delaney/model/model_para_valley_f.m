function p = model_para_valley_f(run)

topopara = 0.05;
p = model_para_valley(topopara);

%%%%%%%%%%%%%%%%%%%
p.options =[-6
            -3
            0
            3
            6];


%% physical
nil = p.coords * 0;  % zeros for all nodes

p.dt = p.options(run); % temperature parameter value... to change with model run.

p.ddf = 0.01/86400; % degree day factor.
p.lr = -0.0075;
p.basal = 7.93e-11 * p.width; % basal melt set to situation A1
p.M = @(t) double(seasonal_source_d(t, p)); % double() is a trick to make anonymous functions faster


%% Initial conditions
p.S_IC = 1 + nil(1:end-1);

%% Time
p.tspan = [0:1:365*3-1]*p.day;
