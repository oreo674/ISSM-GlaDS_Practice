function p = proc_para(p)
% para = proc_para(para)
%
% Calculates the dependent parameters.

% This m-file part of https://bitbucket.org/maurow/1dhydro
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

%% Finite element operators
[p.int_nn, p.int_ne, p.int_ee, p.int_ne_bdy, p.mean_en, p.Dx_en] ...
    = sparseFEM(p.connect, p.coords, p.neumann_nodes);

%% physical parameters
p.phi_m = p.rho_w*p.g*p.B;  % hyd. pot. of bed
p.phi_0 = p.phi_m + p.rho_i*p.g*p.H; % hyd. pot. of ice overburden

%% calculate a good IC for phi: a consistent phi for p.sigma==0
if ~isfield(p, 'phi_IC')
    phi_guess = p.phi_0(p.anodes)/10;
    % this is for sigma==0 as sigma only features in the mass matrix:
    obj = @(phi) double(objectivefun_phi(p.tspan(1), phi, p));

    % some of the options are for Matlab, some for Octave
    opts = optimset('fsolve');
    opts.AutoScaling = 'on';
    opts.TolX = 1e-12;
    opts.TolFun = 1e-12;
    opts.Display = 'off';

    phi = fsolve(obj, phi_guess, opts);
    p.phi_IC = phi;
end
p.IC = [p.S_IC; p.phi_IC];

%% the mass matrix is not diagonal for finite elements:
mass_phi = -bsxfun(@times, p.int_nn(p.anodes,p.anodes), p.sigma(p.anodes));
mass = blkdiag(speye(length(p.S_IC)), mass_phi);

%% ODE
% fixed options:
if p.octave
% $$$     Options for DASPK include:
% $$$
% $$$   keyword                                             value
% $$$   -------                                             -----
% $$$   absolute tolerance                                  1.49012e-08
% $$$   relative tolerance                                  1.49012e-08
% $$$   compute consistent initial condition                0
% $$$   use initial condition heuristics                    0
% $$$   initial condition heuristics
% $$$
% $$$      5.00000
% $$$      6.00000
% $$$      5.00000
% $$$      0.00000
% $$$      0.00000
% $$$      0.01000
% $$$
% $$$   print initial condition info                        0
% $$$   exclude algebraic variables from error test         0
% $$$   algebraic variables                                 0
% $$$   enforce inequality constraints                      0
% $$$   inequality constraint types                         0
% $$$   initial step size                                   -1
% $$$   maximum order                                       5
% $$$   maximum step size                                   -1

    p.mass = mass;
    daspk_options('absolute tolerance', 1.49012e-08);
    daspk_options('relative tolerance', 1.49012e-08);
    daspk_options('enforce inequality constraints', 3);
    daspk_options('inequality constraint types', [p.S_IC*0+1; p.phi_IC*0+0]);

    if any(p.sigma==0)
        p.DAE = true;
        daspk_options('algebraic variables', diag(mass)==0);
        daspk_options('compute consistent initial condition', 1);
    else
        p.DAE = false;
        daspk_options('algebraic variables', 0);
        daspk_options('compute consistent initial condition', 0);
    end
else
    opts = p.odeopts;
    if any(p.sigma==0)
        opts = odeset(opts, 'MassSingular', 'yes');
        p.DAE = true;
    else
        opts = odeset(opts, 'MassSingular', 'no');
        p.DAE = false;
    end
    %% progress updates during runs:
    progress_anon = @(t,y,flag) double(progress(t,y,flag, p.tspan(1), p.tspan(end)));


    %% Jacobian sparsity pattern
    % The sparsity pattern of the Jacobian.  This is optional but speeds up Matlab's finite
    % difference calculation for larger problems.  Better still would be to provide a
    % function which calculates the Jacobian, but that is tedious and error-prone

    % dF/dy consists of these blocks:
    % $$$ ( dF_S/dS   | dF_S/dphi   )
    % $$$ ( dF_phi/dS | dF_phi/dphi )
    nphi = sum(p.anodes);
    nS = length(p.S_IC);

    % S depends only on its local value
    JPattern_S = speye(nS);
    % phi depends on its neighbouring phis and its local value
    B = [ones(nphi,1), ones(nphi,1), ones(nphi,1)];
    JPattern_phi = spdiags(B, -1:1, nphi, nphi);
    % S depends on the gradient of phi:
    B = [ones(nphi+nS,1), ones(nphi+nS,1)];
    if p.BCtype(1)==1
        JPattern_FS_phi = spdiags(B, -1:0, nS, nphi);
    else
        JPattern_FS_phi = spdiags(B, 0:1, nS, nphi);
    end
    % phi depends on S
    B = [ones(nphi+nS,1), ones(nphi+nS,1)];
    if p.BCtype(1)==1
        JPattern_Fphi_S = spdiags(B, 0:1, nphi, nS);
    else
        JPattern_Fphi_S = spdiags(B, -1:0, nphi, nS);
    end
    JPattern = [JPattern_S ,     JPattern_FS_phi;
                JPattern_Fphi_S, JPattern_phi];

    %% set options
    opts = odeset(opts, 'Mass', mass, 'MStateDependence', 'none', ...
                  'OutputFcn', progress_anon, 'JPattern' , JPattern);
    tic; % start ticking
    p.odeopts = opts;
end


end

function status = progress(tt,unused,flag, start_t, end_t)
    % To print status updates during integration for ode15s
    status = 0;
    secs_per_day = 86400;
%     wallt = toc;
%      if length(tt)>0
%          disp([num2str(tt(1)/secs_per_day, '%.1f'), '/', num2str(end_t/secs_per_day, '%.1f'), ' days/days; walltime: ',  num2str(wallt, '%.1f')])
%          tic;
%      end
 end

function dphidt = objectivefun_phi(t, phi, p)
    % objectivefun for phi only; get it out of the objectivefun.m.
    y = zeros(p.n_nodes-1 + length(phi), 1);
    y(1:p.n_nodes-1) = p.S_IC;
    y(p.n_nodes:end) = phi;
    dydt = objectivefun(t, y, p);
    dphidt = dydt(p.n_nodes:end);
end