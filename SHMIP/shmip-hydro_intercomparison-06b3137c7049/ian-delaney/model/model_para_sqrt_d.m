function p = model_para_sqrt_d(run)

p = model_para_sqrt();

%%%%%%%%%%%%%%%%%%%
p.options =[-4 % delta T
            -2
            0
            2
            4];

p.dt = p.options(run); % temperature parameter value... to change with model run.

%% physical
nil = p.coords * 0;  % zeros for all nodes

% otherwise the solver crashes
p.sigma = nil + 5e-12; % Volume storage capacity per unit length per unit pressure [m^2/Pa]

p.ddf = 0.01/86400; % degree day factor.
p.lr = -0.0075;
p.basal = 7.93e-11 * p.width; % basal melt set to situation A1
p.M = @(t) double(seasonal_source_d(t, p)); % double() is a trick to make anonymous functions faster

%% Initial conditions
p.S_IC = 0.1 + nil(1:end-1); %need to change IC so A1 (winter discharge) is met.

%% Time
p.tspan = [0:1:365*3-1]*p.day;
