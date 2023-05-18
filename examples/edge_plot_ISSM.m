function edge_plot_ISSM(ax,md,Cdata,cmap,caxis,varargin)
% edge_plot. Plot variable values that are stored on edges for ISSM-GlaDS.
%
% edge_plot(ax, md, Cdata, cmap, caxis) uses the axes instance ax, the
% model md. Cdata is the plot data used to color edges. cmap is the
% colormap (N x 3 array).
%
% caxis can be 'none' or an array [cmin, cmax]. If 'none', uses simple
% rules to guess good colourbar limits.
%
% edge_plot(ax, md, Cdata, cmap, caxis, 'Lineweights', lw, 'vmin', vmin)
%
% lw is a pair [lwmin, lwmax] giving the minimum and maximum lineweights.
% vmin is a minimum threshold so that only values above vmin are shown.
%
% See also element_plot

% Convert to bamg mesh if necessary
if isnan(md.mesh.edges)
    md = meshconvert(md, md.mesh.elements, md.mesh.x, md.mesh.y);
end

% Default plot parameters
lw_min = 0.5;
lw_max = 2;
% Parse varargin to update default parameters
plot_args = {};
vmin_nan = nan;
for ii=1:2:length(varargin)
   if strcmp(varargin{ii}, 'Lineweights')
       lineweights = varargin{ii+1};
       lw_min = lineweights(1);
       lw_max = lineweights(2);
   elseif strcmp(varargin{ii}, 'vmin')
       vmin_nan = varargin{ii+1};
   end
end

if strcmp(caxis, 'none')
    % Calculate our own automatic bounds on caxis
    cmax = max(Cdata);
    cmin = min(Cdata);

    if cmin>=0
        cmax_scale = 10^(floor(log10(abs(cmax)))-1);
        vmax = cmax_scale*ceil(cmax/cmax_scale);
        if ~isnan(vmin_nan)
            vmin = vmin_nan;
        else
            vmin = 0;
        end
    else
        cmax_scale = 10^(floor(log10( max([abs(cmax), abs(cmin)])))-1);
        vmax = cmax_scale * ceil( max([abs(cmax), abs(cmin)])/cmax_scale);

        if ~isnan(vmin_nan)
            vmin = vmin_nan;
        else
            vmin = - vmax;
        end
    end
else
    vmin = caxis(1);
    vmax = caxis(2);
end
getix=@(v) min( max(1, round(size(cmap,1)*(v-vmin)/(vmax-vmin) )), size(cmap,1));
getcolour=@(v) cmap(getix(v),:);

getlw=@(v) min( max(lw_min, (lw_max*(v-vmin)/(vmax-vmin) )), lw_max);

lws = getlw(Cdata);
linecolours = getcolour(Cdata);

hold on
for ii=1:md.mesh.numberofedges
    nodesx = md.mesh.x(md.mesh.edges(ii, 1:2));
    nodesy = md.mesh.y(md.mesh.edges(ii, 1:2));
    if isnan(Cdata(ii))
        linecolours(ii) = [0,0,0,0];
        lws(ii) = 0.01;
    end
    line(ax,nodesx,nodesy,'color',linecolours(ii,:),'linewidth',lws(ii));
end

colormap(ax,cmap);
cb=colorbar(ax);
ticks=linspace(0,1,5);
tickl=linspace(vmin,vmax,5);
ticklabels=cellstr(num2str(tickl'));
set(cb,'Ticks',ticks,'TickLabels',ticklabels)
axis image
