function [Q, Xi, creep, melt, N, grad_phi, channelized, stored, tot_stored, source, tot_source, tau, dm_max, x, x2] =postproc(phi, S, p)
% [Q, Xi, creep, melt, N, grad_phi, channelized, stored, tot_stored, source, tot_source, x, x2] =postproc(phi, S, para);
%
% Calculates a variety of interesting quantities from the model solution.
%
% Q  -- discharge [m^3/s]
% Xi -- stream power [W/m]
% creep -- creep closure rate [m^2/s]
% melt -- melt opening term (note, there so no pressure-melt term) [m^2/s]
% N -- effective pressure [Pa]
% grad_phi -- gradient of phi [pa/m]
% channelized -- Schoof 2010 criterion (NOT WORKING)
% stored -- stored water per unit channel length [m^2]
% tot_stored -- total stored water [m^3]
% source -- source term per unit channel length [m^2/s]
% tot_source -- total source [m^3/s]
% tau -- shear stress (?)
% dm_max -- maximal transportable bedload size
% x -- node coordinates [m]
% x2 -- cell center coordinates [m]

% This m-file part of https://bitbucket.org/maurow/1dhydro
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

x = p.coords;
x2 = (x(1:end-1)+x(2:end))/2;

for t=1:size(phi,2)
    grad_phi(:,t) = p.Dx_en * phi(:,t);
    Q(:,t)  = - p.k * S(:,t).^p.alpha .* (grad_phi(:,t).^2+p.tiny^2).^((p.beta-2)/2) .* grad_phi(:,t) ;
    Xi(:,t) = - grad_phi(:,t).*Q(:,t);
    N(:,t)  = p.mean_en * (p.phi_0 - phi(:,t));  % N on elements
    creep(:,t) = p.A*S(:,t).*N(:,t).^p.n;
    melt(:,t) = -Xi(:,t)/p.L/p.rho_i;
    % critical discharge, see Schoof 2010
    % THIS IS NOT WORKING AS I WOULD EXPECT?!
    c1 = 1/(p.L*p.rho_i);
    channelized(:,t) = (p.u*p.h ./ (c1*(p.alpha-1)*grad_phi(:,t)) ) < Q(:,t);
    % storage
    pw = phi(:,t)-p.phi_m;
    stored(:,t) = p.sigma.*pw; % [m^2]
    tot_stored(t) = sum(p.sigma.*pw.*gradient(x));
    if isa(p.M, 'function_handle')
        source(:,t) = p.M(p.tspan(t));
    else
        source(:,t) =  p.M;
    end
    tot_source(t) = sum(source(:,t).*gradient(x));

end
tau=1/8*0.12*p.rho_w*(Q./S).^2;
dm_max = tau./(p.mu .* p.g .*(p.rho_s-p.rho_w));
