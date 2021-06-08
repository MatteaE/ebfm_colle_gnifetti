%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C] = func_init_constants()

%% Tuning parameters (to be calibrated against observational data if available)
C.alb_fresh = 0.83;         % albedo fresh snow
C.alb_firn = 0.52;          % albedo firn
C.alb_ice = 0.39;           % albedo ice
C.tstar_wet = 10;           % albedo decay time-scale (wet snow)
C.tstar_dry = 30;           % albedo decay time-scale (dry snow)
C.tstar_K = 14;              % albedo decay time-scale coefficient
C.dstar = 7.0;              % albedo characteristic depth (mm w.e.)
C.rainsnowT = 273.75;       % temperature of snow to rain transition (K)
C.b = 0.420;                % constant in LWin formulation
C.ecl = 0.960;              % constant in LWin formulation
C.Dfreshsnow = 350.0;       % density of fresh snow (kg m-3)
C.lambda = 1.6;             % optical thickness empirical constant
C.k_aer = 0.97;             % aerosol transmissivity exponent
C.turb = 0.0025;            % background turbulent exchange coefficient
C.Pthres = 2.5d-8;          % threshold precipitation to reset time since last snowfall (m w.e. s-1)
C.Trunoff = 0.001;          % slush runoff time-scale (days)
C.geothermal_flux = 0.04;   % geothermal heat flux (W m-2)
C.soil_albedo = 0.15;       % soil albedo

% C.Ch = 0.00127;             % turbulent exchange coefficient for bulk approach (Oerlemans 2000).
C.dz = 2.0;                 % vertical extent of the layer used for bulk turbulent fluxes.
C.z0 = 0.001;                 % Roughness length for the Essery (2004) approach.
% C.z0u = 0.001;              % surface roughness length for wind speed.
% C.z0t = 0.0001;             % surface roughness length for temperature.
% C.z0q = 0.0001;             % surface roughness length for humidity.


%% Here we compute the (constant) log mean heights used in Oke (1987) formulation of turbulent fluxes.
%C.zu = (C.dz - C.z0u) / (log(C.dz / C.z0u));
%C.zt = (C.dz - C.z0t) / (log(C.dz / C.z0t));
%C.zq = (C.dz - C.z0q) / (log(C.dz / C.z0q));

%% Essery bulk approach.
C.Chn = 0.16 * (log(C.dz / C.z0)^-2);         % neutral exchange coefficient for Essery (2004) bulk approach.
C.fz = 0.25 * sqrt(C.z0 / C.dz);              % parameter in the Essery (2004) bulk approach.


%% Physical constants
C.p = 2;                    % exponent in LWin formulation
C.boltz = 5.67d-8;          % Stefan-Boltzmann constant (W m-2 K-4)
C.VP0 = 610.5;              % vapor pressure at 0 degrees C
C.Cp = 1005.7;              % specific heat of dry air (J kg-1 K-1)
C.Cw = 4187.0;              % specific heat of water (J kg-1 K-1)
C.Ls = 2.83d6;              % latent heat of sublimation/riming (J kg-1)
C.Lm = 0.33d6;              % latent heat of melting/fusion (J kg-1)
C.Lv = 2.5d6;               % latent heat of evaporation/condensation (J kg-1)
C.Rv = 462.0;               % specific gas constant water vapor (J kg-1 K-1)
C.Rd = 287.0;               % specific gas constant dry air (J kg-1 K-1)
C.eps = 0.622;              % C.Rd/C.Rv
C.dTacc = 0.01;             % threshold dT in solving the energy balance equation (K)
C.Pref = 1015d2;            % reference air pressure (Pa)
C.Dice = 900.0;             % density of ice (kg m-3)
C.Dwater = 1000.0;          % density of water (kg m-3)
C.g = 9.81;                 % gravitational acceleration (m s-2)
C.Ec = 60000;               % gravitational densification factor
C.Eg = 42400;               % gravitational densification factor
C.rd = 8.314;               % universal gas constant (J mol-1 K-1)          
C.k_turb = 0.0004;          % turbulent flux coefficient
C.Pr = 5;                   % Prandtl number in SHF/LHF formulation
C.T0 = 273.15;              % melting temperature ice (K)
C.soil_Kfrozen = 2.50;      % soil conductivity (frozen)
C.soil_Kthawed = 1.45;      % soil conductivity (thawed)
C.soil_Cfrozen = 1.95d6;    % soil specific heat capacity (frozen)
C.soil_Cthawed = 2.45d6;    % soil specific heat capacity (thawed)
C.soil_THwmin = 0.05;       % soil minimum water content
C.soil_THwmax = 0.30;       % soil maximum water content
C.soil_delta = 0.17;        % soil freezing parameter (deg C)
C.soil_density = 1000.0;    % soil density (no functionality)

C.k = 0.40;                 % von Karman constant

end
