% This runs the model.
% Uncomment different model_para* calls for different scenarios.

% This m-file part of https://bitbucket.org/maurow/1dhydro
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

suites = {'A', 'B', 'C', 'D', 'E', 'F'};
runs = [6,5,4,5,5,5];

for i = 1:length(suites)
    suite = suites{i};
    for run = 1:runs(i)

        %% choose one of the parameter files:
        tic
        switch suite
          % case 'A'
          %   para = model_para_sqrt_a(run);
          % case 'B'
          %   para = model_para_sqrt_b(run);
          % case 'C'
          %   para = model_para_sqrt_c(run);
          case 'D'
            para = model_para_sqrt_d(run);
          % case 'E'
          %   para = model_para_valley_e(run);
          % case 'F'
          %   para = model_para_valley_f(run);
        end
        %% declare run
        run_name = [suite, num2str(run), '_idel']

        para = proc_para(para);

        %% run it

        % The objective function, this is where all the physics is,
        % except some is also in the mass-matrix which is set in proc_para.
        objfun = @(t, y) double(objectivefun(t, y, para));

        %% solve
        [t, out] = ode15s(objfun, para.tspan, para.IC, para.odeopts);
        if t(end)~=para.tspan(end)
            warning(['Run not completed: ', suite, num2str(run), '. Skipping!'])
            delete(['../results/',run_name,'.mat'])
            continue
        else
            para.tspan = t;
        end

        % unwrap
        S = out(:,1:para.n_nodes-1)';
        phi = zeros(para.n_nodes,size(out,1));
        phi(para.anodes,:) = out(:,para.n_nodes:end)'; % on nodes
        phi(~para.anodes,:) = repmat(para.BCval(para.BCtype==1), 1, size(out,1)); % set Dirichlet BC

        %% post-processing
        [Q, Xi, creep, melt, N, grad_phi, channelized, stored, tot_stored,...
         source, tot_source, tau, dm_max, x, x2] = postproc(phi, S, para);
        run_time = toc

        save(['../results/',run_name])
    end
end
