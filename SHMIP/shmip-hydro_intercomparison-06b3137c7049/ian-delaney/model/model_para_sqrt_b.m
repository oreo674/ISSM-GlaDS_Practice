function p = model_para_sqrt_b(run)

p = model_para_sqrt();

x1 = p.coords(end);
dx = p.coords(2)-p.coords(1);
num = length(p.coords);

%% physical
nil = p.coords * 0;  % zeros for all nodes

mou = dlmread(['../../input_functions/source/B',num2str(run),'_M.csv'],',');

p.M = p.coords*0 ; %basal
for i= 1:size(mou,1) % has trouble with B1...
    for j = 1:length(p.M)
        if mou(i,2) == j*x1/num  %find node where moulin is.
            p.M(j) = mou(i,4)/dx+p.M(j);
        end
    end
end
if abs(sum(p.M)-sum(mou(:,4)/dx))>1e-5
    error('Not all moulins included')
end

% add basal melt:
p.M = p.M + (7.93e-11 * p.width) ; % is source term positive or negative

% add a small sigma to help the ode15s solver
p.sigma = nil + 1e-12;


%% Initial conditions
p.S_IC = 0.5 + nil(1:end-1);

%% Time
p.tspan = [0:10:500]*p.day;
