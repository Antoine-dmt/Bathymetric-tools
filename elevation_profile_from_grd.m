clear all
%script for loading a .grd file and measure elevation lines on it
%using mapping toolbox from matlab
%A. Demont & J.A. Olive & J. Escartin 27/04/21

[X,Y,Z] = grdread2('SATW.grd');

precision = 1000;

figure
K = surf(X,Y,Z);
set(K,'Linestyle','none')
view(2)
profile_boundaries = ginput;

%profile_boundaries = [-45.4288 23.5135; -45.1584 23.4691];

vect_lon = linspace(profile_boundaries(1,1), profile_boundaries(2,1),precision);
vect_lat = linspace(profile_boundaries(1,2), profile_boundaries(2,2),precision);


%interpolating data
Zp = interp2(X,Y,Z,vect_lon,vect_lat);
alen = distance(profile_boundaries(1,2),profile_boundaries(1,1),profile_boundaries(2,2),profile_boundaries(2,1));
d = deg2km(alen);
dist = linspace(0,d, precision)*1e3;

%profile plot
figure
plot(dist,Zp,'.')
%axis equal
title('Elevation profile from map')
ylabel('Depth (m)')
xlabel('Distance (m)')

%%%Calculating OCC radius and fault shallow dip
array_positions = ginput;% 6 positions expected, break away top, break away bottom,OCC top, OCC bottom, fault top dip sample, fault bottom dip sample 

top_break_away = array_positions(1,2);
bot_break_away = array_positions(2,2);
OCC_top = array_positions(3,2);
OCC_bot = array_positions(4,2);
fault_top = array_positions(5,1);
fault_bot = array_positions(6,1);

break_height = top_break_away - bot_break_away;
height_OCC = abs(OCC_top - OCC_bot);


%%%%Curve radius OCC routine

fit_dist = dist((dist > array_positions(3,1) & dist< array_positions(4,1))); %fitting topography
fit_elev = Zp(dist > array_positions(3,1) & dist< array_positions(4,1));

%fitting now curvature radius
fit_deg = 2; %degree for polynomial fitting
polyfit_curve  = polyfit(fit_dist ,fit_elev, fit_deg);
x1 = fit_dist(1):1:fit_dist(end);
y1 = polyval(polyfit_curve,x1);
figure
plot(fit_dist,fit_elev,'o',x1,y1,'--')


%grad_first  = gradient(y1(:)) ./ gradient(x1(:));
%curvature_array = gradient(grad_first(:)) ./ gradient(x1(:));
%radius_array = 1./ (abs(curvature_array));
%radius_array(isinf(radius_array)) = NaN;
radius_average =abs( 1/ (2*polyfit_curve(1)));

%finding dip at fault emergence
%not true fault but fault sampling by OCC shape considered as a valid true
%fault dip proxy
Xfault = dist(dist>= fault_top & dist <= fault_bot);
Yfault = Zp(dist>= fault_top & dist <= fault_bot); 

polyfit_curve_dip = polyfit(Xfault,Yfault,2);
%yfit_fault = polyval(polyfit_curve
x_shallow = max(Xfault);
dip_shallow_value = atand(abs(x_shallow*polyfit_curve_dip(1)*2 +polyfit_curve_dip(2)));

save('SATW_morpho_data','dist','Zp','break_height','height_OCC','radius_average','dip_shallow_value')
