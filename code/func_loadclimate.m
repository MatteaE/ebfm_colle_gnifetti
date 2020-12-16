%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load climate fields from file or heuristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [clim,A] = func_loadclimate(C,grid,clim,io,t,time,A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read climate input from file
% - in this example time-dependent, but spatially fixed values are read for 
%   temperature, relative humidity, cloud cover, air pressure and
%   precipitation
% - elevation corrections are applied for temperature, precipitation and
%   air pressure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (t==1)
    load(io.climfile);
    A.data_clim = climate;
end

tt = round((time.TCUR - time.TI)/time.dt) + 1;

%% Reference elevation for climate fields [m a.s.l.]
zref = 4560; 

%% Temperature [K]
% - Decreases with elevation using a lapse rate
clim.T(:) = A.data_clim(tt,7) + A.data_clim(tt,8)*(grid.z_mask-zref);

%% Precipitation [m w.e. timestep-1]
% - Increases with elevation (see Van Pelt et al. 2019)
% P_0 = 1.11;
% P_1 = 0.0022;
% clim.P(:) = A.data_clim(tt,9) .* max(P_0 + (min(grid.z_mask,900)-zref) * P_1,0.1);
clim.P(:) = A.data_clim(tt,9) .* grid.accum_clim  .* A.data_clim(tt,15);

%% Cloud cover [fraction]
clim.C(:) = A.data_clim(tt,10);

%% Relative humidity [fraction]
clim.RH(:) = A.data_clim(tt,11);

%% Air pressure [Pa]
clim.Pres(:) = A.data_clim(tt,12) * exp(A.data_clim(tt,13)*(grid.z_mask-zref));

%% Potential temperature lapse rate [K m-1]
clim.Theta_lapse = A.data_clim(tt,14);

%% Wind speed [m s-1]
clim.Wind_speed = A.data_clim(tt,16);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Derived climate fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Snowfall / Rainfall
clim.snow = clim.P .* (clim.T < C.rainsnowT-1);
clim.rain = clim.P .* (clim.T > C.rainsnowT+1);
clim.snow = clim.snow + clim.P .* (C.rainsnowT-clim.T+1)./2 .* (clim.T < C.rainsnowT+1 & clim.T > C.rainsnowT-1);
clim.rain = clim.rain + clim.P .* (1+clim.T-C.rainsnowT)./2 .* (clim.T < C.rainsnowT+1 & clim.T > C.rainsnowT-1);

%% Annual snow accumulation (needed for gravitational settling calculations)
A.ys = (1.0-1.0/(365.0/time.dt)).*A.ys + clim.P.*1d3;
logys = log(A.ys);
clim.yearsnow = repmat(A.ys,[1 grid.nl]);
clim.logyearsnow = repmat(logys,[1 grid.nl]);

%% Vapor pressure (VP) / specific humidity (q)
VPsat = C.VP0.*exp(C.Lv/C.Rv.*(1.0./273.15-1.0./clim.T)) .* (clim.T>=273.15) + ...
        C.VP0.*exp(C.Ls/C.Rv.*(1.0./273.15-1.0./clim.T)) .* (clim.T<273.15);
clim.VP = clim.RH .* VPsat;
clim.q = clim.RH .* (VPsat .* C.eps ./ clim.Pres);

%% Saturation vapor pressure and humidity at the snow surface (Essery and Etchevers 2004).
clim.VPsat_surf = C.VP0.*exp(C.Lv/C.Rv.*(1.0./273.15-1.0./A.Tsurf)) .* (A.Tsurf>=273.15) + ...
C.VP0.*exp(C.Ls/C.Rv.*(1.0./273.15-1.0./A.Tsurf)) .* (A.Tsurf<273.15);
clim.qsat_surf = (clim.VPsat_surf .* C.eps ./ clim.Pres);
    
%% Air density (Dair)
clim.Dair = clim.Pres./C.Rd./clim.T;

%% Time since last snow fall event (timelastsnow)
A.timelastsnow(clim.snow/(time.dt*24*3600)>C.Pthres) = time.TCUR;

%% Potential temperature (Theta)
clim.Theta = clim.T.*(C.Pref./clim.Pres).^(C.Rd/C.Cp);



A.climT = clim.T;
A.climP = clim.P;
A.climC = clim.C;
A.climRH = clim.RH;
A.climPres = clim.Pres;
A.climsnow = clim.snow;
A.climrain = clim.rain;
A.climWindspeed = clim.Wind_speed;

%A.climRi = clim.Ri;
A.climRi2 = clim.Ri2;


end
