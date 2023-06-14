function [bed_ele, ice_thickness] = sqrttopo(x, y, bedslope)
% [bed_ele, ice_thickness] = sqrttopo(x, y, bedslope)
%
% Surface topography with square root shape and linearly sloping bed.
%
% Input:
% - slope: positive means increasing bed elevation

bed_ele = x*bedslope;

thick_para = 6;
thick_off = 5e3;
thick_nz = 1;
surf = thick_para*sqrt(x+thick_off) - thick_para*sqrt(thick_off) + thick_nz;
ice_thickness = surf-bed_ele;
