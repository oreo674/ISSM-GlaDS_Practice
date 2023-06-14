function [int_nn, int_ne, int_ee, int_ne_bdy, mean_en, Dx_en] = sparseFEM(connect, x, neumann_nodes);
% [int_nn, int_ne, int_ee, int_ne_bdy, mean_en, Dx_en] = sparseFEM(connect, coords, neumann_nodes);
%
% Sets up integration, averaging and differentiation operators for piecewise linear finite
% elements in 1 dimensions
%
% Example usage for residual calculation of Poisson equation with finite elements:
%     res =  (Dx_en' *int_ee* Dx_en) * phi ...
%           - int_nn*f - int_ne_bdy * neumann;
% 
% INPUT
%   connect:          n_elements-by-(dimension+1) array of vertex indices for each element
%   coords:           dimension-by-n_nodes array of Cartesian coordinates of nodes
%   neumann_nodes:    index of Neumann BC nodes
%
% OUTPUT
%
% n --> node
% e --> element
%
% Normal element operators:
%   int_nn: f.'*int_nn*g integrates the product of two functions, both defined by linear
%           interpolation between nodal values f and g
%   int_ne: f.'*int_ne*h integrates the product of two functions, one defined by
%           interpolation between nodal values f, the other by piecewise constant
%           values h on elements
%   int_ee: k.'*int_ee*h integrates the product of two functions, each piecewise
%           constant values k and h on elements
%   mean_en:    mean_en*f computes mean over each element of a function defined by
%               linear interpolation between nodal values f
%   Dx_en:  Dx_en*f computes on each element the x-derivative of  a function defined by
%           linear interpolation between nodal values f
%   int_nf_bdy: f.'*int_nf_bdy*h integrates over the boundary the product of two functions, one defined by
%               interpolation between nodal values f, the other by piecewise constant
%               values h on boundary factes. (here in 1D this integral is trivial as the facet is a point!)

% This m-file is modified from the sparseFEM package https://bitbucket.org/maurow/sparsefem
% Copyright (c) 2014, Christian Schoof, Ian J Hewitt & Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

n_nodes = size(x,1);                 %number of nodes
n_elements = size(connect,1);           %number of elements

elements = (1:n_elements)';
facets_bdy = (1:length(neumann_nodes))';

%operators defined on entire mesh
connect1 = connect(:,1);
connect2 = connect(:,2);
dx_connect = x(connect2)-x(connect1); % Jacobian
abs_dx_connect = abs(dx_connect);   % element size
%exact quadrature
int_nn = sparse([connect1;           connect1;           connect2;           connect2], ...
                [connect1;           connect2;           connect1;           connect2], ...
                [1/3*abs_dx_connect; 1/6*abs_dx_connect; 1/6*abs_dx_connect; 1/3*abs_dx_connect], n_nodes,n_nodes);

int_ne =  sparse([connect1; connect2], [elements; elements], [1/2*abs_dx_connect; 1/2*abs_dx_connect], n_nodes,n_elements);
int_ee =  sparse(elements, elements, abs_dx_connect, n_elements, n_elements);
mean_en = sparse([elements; elements], [connect1; connect2], 1/2*ones(2*n_elements,1), n_elements,n_nodes);     
Dx_en =   sparse([elements; elements], [connect1; connect2], [-1./dx_connect; 1./dx_connect], n_elements,n_nodes);
    
int_ne_bdy  = sparse(neumann_nodes, facets_bdy, ones(size(neumann_nodes)), n_nodes, length(neumann_nodes));
