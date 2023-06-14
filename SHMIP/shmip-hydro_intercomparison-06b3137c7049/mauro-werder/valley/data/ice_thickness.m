function out = ice_thickness(xy, time, bed_para)

[~,out] = valley(xy(:,1), xy(:,2), bed_para);
