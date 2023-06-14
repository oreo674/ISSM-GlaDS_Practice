%% Script to save set of model runs to .nc files
% files from all suites and runs should be in a directory (../../results)
% to run correctly. .nc files are writen to a directory
% (../../results/netCDF_files/)

suites = {'A', 'B', 'C', 'D', 'E', 'F'};
runs = [6,5,4,5,5,5];

for i = 1:length(suites)
    suite = suites{i};
    for j = 1:runs(i)

        load(['../results/',suite, num2str(j),'_idel.mat'],'-mat');

        filename = ['../results/netCDF_files/',run_name,'.nc']; % make file name compatible with naming convension
        disp(run_name)
        % A series- steady
        switch suite
          case {'A', 'B', 'E'}
            inds = length(para.tspan);
          case 'C'
            inds = length(para.tspan)-23:length(para.tspan);
          case {'D', 'F'}
            inds = length(para.tspan)-364:length(para.tspan);
        end
        time = para.tspan(inds) - para.tspan(inds(1));
        ntime = length(time);

        % NETCDF dimensions
        dim = 1;
        index_ch= length(x2);
        index1=length(x);

        n_nodes_ch=2;
        dims = { 'time', length(time); % / time dimension
                 'dim', dim;% / spatial dimensions
                 'n_nodes_ch', n_nodes_ch ; % nodes per channel
                 'index1', index1; % / index into coords1
                 'index_ch', index_ch; % / channel index
               };


        mySchema.Name   = '/';
        mySchema.Format = 'classic';
        for i = 1:size(dims,1)
            n = dims{i,1};
            v = dims{i,2};
            mySchema.Dimensions(i).Name   = n;
            mySchema.Dimensions(i).Length = v;
        end
        ncwriteschema(filename, mySchema);

        %% create dimensional variable time
        nccreate(filename, 'time', 'Dimensions',{'time' ntime});
        ncwrite(filename, 'time', time);
        ncwriteatt(filename, '/time', 'units', 's')
        ncwriteatt(filename, '/time', 'long_name', 'time')

        %% create and write physical variables

        % create and write coords1
        nccreate(filename, 'coords1', 'Dimensions',{'index1' index1, 'dim' dim });  % might need to flip
        ncwrite(filename, 'coords1', x);
        ncwriteatt(filename,'/coords1', 'long_name', 'node coordinates');
        ncwriteatt(filename,'/coords1', 'units', 'm');

        % create and write coords_ch
        nccreate(filename, 'coords_ch', 'Dimensions',{'index_ch' index_ch, 'dim' dim });  % might need to flip
        ncwrite(filename, 'coords_ch', x2);
        ncwriteatt(filename,'/coords_ch', 'long_name', 'channel midpoint coordinates');
        ncwriteatt(filename,'/coords_ch', 'units', 'm');

        % create and write connectivity of channels
        nccreate(filename, 'connect_ch', 'Dimensions',{'index_ch' index_ch, 'n_nodes_ch' n_nodes_ch});  % might need to flip
        ncwrite(filename, 'connect_ch',[0:length(x)-2;1:length(x)-1]');
        ncwriteatt(filename,'/connect_ch', 'long_name', 'channel connectivity');
        ncwriteatt(filename,'/connect_ch', 'units', '');

        % create and write ice thickness, H
        nccreate(filename, 'H', 'Dimensions',{ 'index1'});
        ncwrite(filename, 'H', para.H);
        ncwriteatt(filename,'/H', 'long_name', 'ice thickness');
        ncwriteatt(filename,'/H', 'units', 'm');

        % create and write bed elevation, B
        nccreate(filename, 'B', 'Dimensions',{'index1'});
        ncwrite(filename, 'B', para.B);
        ncwriteatt(filename,'/B', 'long_name', 'bed elevation');
        ncwriteatt(filename,'/B', 'units', 'm')

        % create and write bed elevation, W
        nccreate(filename, 'W', 'Dimensions',{'index1'});
        if length(para.width) == 1
            ncwrite(filename, 'W', repmat(para.width, 1, length(x))); % there is a single value for para.width in sqrt glaciers...
        else
            ncwrite(filename, 'W', para.width);
        end
        ncwriteatt(filename,'/W', 'long_name', 'domain width');
        ncwriteatt(filename,'/W', 'units', 'm')

        % create and write effective pressure, N
        nccreate(filename, 'N', 'Dimensions',{'index1' index1, 'time' ntime });
        ncwrite( filename, 'N', N(:,inds));
        ncwriteatt(filename,'/N', 'long_name', 'effective pressure');
        ncwriteatt(filename,'/N', 'units', 'Pa')

        % create and write stored water effective layer thickness, hstore
        nccreate(filename, 'hstore', 'Dimensions',{'index1' index1, 'time' ntime});
        ncwrite(filename, 'hstore',  bsxfun(@rdivide, stored(:,inds), para.width));
        ncwriteatt(filename,'/hstore', 'long_name', 'stored water effective layer thickness');
        ncwriteatt(filename,'/hstore', 'units', 'm')

        % create and write channel cross-section, S
        nccreate(filename, 'S', 'Dimensions',{ 'index_ch' index_ch,'time' ntime});
        ncwrite(filename, 'S', S(:,inds));
        ncwriteatt(filename,'/S', 'long_name', 'channel cross-sectional area');
        ncwriteatt(filename,'/S', 'units', 'm^2')

        % create and write channel discharge, Qs
        nccreate(filename, 'Q', 'Dimensions',{'index_ch' index_ch, 'time' ntime});
        ncwrite(filename, 'Q', abs(Q(:,inds)));
        ncwriteatt(filename,'/Q', 'long_name', 'channel discharge');
        ncwriteatt(filename,'/Q', 'units', 'm^3/s')

        %% global attributes
        ncwriteatt(filename, '/', 'title', ['delaney_',run_name])
        ncwriteatt(filename, '/', 'meshtype', 'structured')
        ncwriteatt(filename, '/', 'channels_on_edges', 'yes')
        ncwriteatt(filename, '/', 'institution', 'Ian Delaney, ETH-Zuerich')
        ncwriteatt(filename, '/', 'source', 'https://bitbucket.org/maurow/1dhydro/overview')
        ncwriteatt(filename, '/', 'references', 'http://shmip.bitbucket.org/')

        %ncdisp(filename)
    end
end