function dydt = objectivefun(t, y, p)
% dydt = objectivefun(t, y, para);
%
% Calculates dS_dt and dphi_dt for use with ode15s.

% This m-file part of https://bitbucket.org/maurow/1dhydro
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

%% setup
% split y into its components phi and S:
S = y(1:p.n_nodes-1);     % on elements
phi = zeros(p.n_nodes,1);
phi(p.anodes) = y(p.n_nodes:end); % on active nodes
phi(~p.anodes) = p.BCval(p.BCtype==1); % set Dirichlet BC

if any(S<0)
    error('S smaller than zero!')
end

%% update source term if it is a function
if isa(p.M, 'function_handle')
    M = p.M(t);
else
    M = p.M;
end

%% calculate dphi_dt (or residual of phi equation if there is no storage).
% This uses a Finite Element discretisation:
grad_phi = p.Dx_en * phi;
Q  = - p.k * S.^p.alpha .* (grad_phi.^2+p.tiny^2).^((p.beta-2)/2) .* grad_phi ;
Xi = - grad_phi.*Q;
N  = p.mean_en * (p.phi_0 - phi);  % N on elements

% calc dphi_dt
dphi_dt = -p.Dx_en(:,p.anodes)' *p.int_ee* Q + ...
           p.int_ne(p.anodes,:) * (Xi/p.L * (1/p.rho_i - 1/p.rho_w) ...
                                   + p.u*p.h - p.A*S.*N.^p.n - p.mean_en*M);
% set Neumann BC
dphi_dt = dphi_dt + p.int_ne_bdy(p.anodes,:) * p.BCval(p.BCtype==2);

%% calculate dS/dt
dS_dt = Xi/(p.rho_i*p.L) + p.u*p.h - p.A*S.*N.^p.n;
% Note, no boundary conditions are needed as this is an ODE

%% put dS_dt and dphi_dt into dydt
dydt = [dS_dt; dphi_dt];

if any(imag(dydt))
    error('Imaginary values!')
end
