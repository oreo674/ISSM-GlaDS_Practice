function p = model_para_sqrt_a(run)

p = model_para_sqrt();

p.options = [7.93e-11 % input options of runs: source in [m/s]
             1.59e-9
             5.79e-9
             2.50e-8
             4.50e-8
             5.79e-7];

p.source = p.options(run); % source give various scenarios [m/s]

%% physical
nil = p.coords * 0;  % zeros for all nodes

p.M = p.coords*0 + (p.source * p.width) ; % is source term positive or negative

if run==6
    %% Initial conditions
    p.S_IC = 5 + nil(1:end-1);
end

%% Time
p.tspan = [0:10:500]*p.day;
