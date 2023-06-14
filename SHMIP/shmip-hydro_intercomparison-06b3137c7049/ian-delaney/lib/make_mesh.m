function [connect, coords, n_nodes, anodes, neumann_nodes] = make_mesh(a, b, n, BCtype, rand_fac)
%  [connect, coords] = make_mesh(a, b, n, BCtype, rand_fac)
%
% Makes a 1D mesh
%
% 1 argument call: argument are the vertex positions
%
% 2 or 3 argument call: 
%  - a: start of interval
%  - b: end
%  - n: number of points (100 default)
%  - BCtype: boundary type, 1==Dirichlet, 2==Neumann; default [1,1];
%  - rand_fac: adds a random offset to each node coordinate

% This m-file is modified from the sparseFEM package https://bitbucket.org/maurow/sparsefem
% Copyright (c) 2014, Christian Schoof, Ian J Hewitt & Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

if nargin==1
    coords = a;
    n = length(coords);
else
    if ~exist('a', 'var') || isempty(a)
        a=0;
    end
    if ~exist('b', 'var') || isempty(b)
        b=1;
    end
    if ~exist('n', 'var') || isempty(n)
        n=100;
    end
    coords = linspace(a,b,n)';
end

if ~exist('BCtype', 'var') || isempty(BCtype)
    BCtype = [1,1]; 
end

connect = [(1:n-1)', (2:n)'];

% randomize if requested
if exist('rand_fac', 'var') && ~isempty(rand_fac)
    if rand_fac>=0.5
        error('rand_fac needs be smaller than 0.5')
    end
    dx = coords(2)-coords(1);
    coords = coords + rand(n,1)*rand_fac*dx;   
end

%% BC
n_nodes = length(coords);
anodes = logical(ones(size(coords))); 
anodes([1,end]) = BCtype==2;          % only Neumann boundary nodes are active

tmp = [1, n_nodes]';
neumann_nodes = tmp(BCtype==2);
