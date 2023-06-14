% This output the tuning data as ascii file
% N, q, Q

meshnr = 4;

% output locations
binlen = 2500;
x = binlen/2:binlen:100e3;
bins = 0:binlen:100e3;
nbins = length(bins)-1;

for scenario=[3,5]
    steady_out = load(['outputmat/sqrt_ch',num2str(scenario),'_mesh',num2str(meshnr),'.mat']);
    [dis, vel, dis_abs, vel_abs, elements, connect_el, dis_c, vel_c, N_eff, N_eff_mean, N_eff_std, N_eff_mean_ts, melt_s, melt_c, melt_tot_ts, mesh, tot_melt_s] = pp_calc_stuff(steady_out);

    qq = dis(:,:,end);
    qq = squeeze(sqrt(qq(:,1).^2+qq(:,2).^2));
    QQ = dis_c;
    %    QQ(QQ<0.01) = NaN;

    N = zeros(nbins,3);
    q = zeros(nbins,3);
    Q = zeros(nbins,1);

    cell_midpoints = zeros(size(steady_out.para.dmesh.tri.connect,1),2);
    for i=1:size(steady_out.para.dmesh.tri.connect,1)
        cell_midpoints(i,:) = sum(steady_out.para.dmesh.tri.nodes(steady_out.para.dmesh.tri.connect(i,:),:))/3;
    end

    for n=1:nbins
        inds = steady_out.para.dmesh.tri.nodes(:,1)>=bins(n) & steady_out.para.dmesh.tri.nodes(:,1)<bins(n+1);
        inds_q = cell_midpoints(:,1)>=bins(n) & cell_midpoints(:,1)<bins(n+1);
        inds_Q = steady_out.para.dmesh.tri.edge_midpoints(:,1)>=bins(n) & steady_out.para.dmesh.tri.edge_midpoints(:,1)<bins(n+1);
        N(n,:) = [mean(N_eff(inds,end)), min(N_eff(inds,end)), max(N_eff(inds,end))];
        q(n,:) = [mean(qq(inds_q,end)), min(qq(inds_q,end)), max(qq(inds_q,end))];
        Q(n) = max(abs(QQ(inds_Q,end)));
    end

    figure
    ax1=subplot(3,1,1)
    plot(x/1e3, N/1e6, 'linewidth',2)
    ylabel('N (MPa)')
    ax2=subplot(3,1,2)
    plot(x/1e3, q, 'linewidth',2)
    ylabel('q (m^2/s)')
    ax3=subplot(3,1,3)
    plot(x/1e3, Q, 'linewidth',2)
    ylabel('Q (m^3/s)')
    xlabel('x (km)')
    linkaxes([ax1,ax2,ax3], 'x')
    axis tight
    print(gcf, '-dpng', ['tuning_A', num2str(scenario), '.png'])

    out = [x',N,q,Q];
    csvwrite(['tuning_A', num2str(scenario)], out)
    % header
    % x,N_mean,N_min,N_max,q_mean,q_min,q_max,Q_max
end
