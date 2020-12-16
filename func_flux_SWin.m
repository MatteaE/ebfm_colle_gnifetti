%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute incoming shortwave radiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SWin,insol] = func_flux_SWin(C,insol,clim)

%% Direct, diffuse and total radiation after shading
insol.TOAdir = (0.2+0.65.*(1-clim.C)).*(1-insol.shade).*insol.TOA;
insol.TOAdiff = (0.8-0.65.*(1-clim.C)).*insol.TOA;
insol.TOAshade = insol.TOAdir + insol.TOAdiff;

%% Transmissivity after gaseous absorption / scattering
m = 35.*(clim.Pres./C.Pref).*(1224.*(insol.TOAflat./insol.I0).^2+1).^(-0.5);
t_rg = 1.021-0.084.*sqrt(m.*(949d0*(clim.Pres./1d3).*1d-5+0.051));

%% Transmissivity after water vapor absorption
tempdew_kelvin = (1/273.15-(C.Rv/C.Ls).*log(clim.q.*clim.Pres./(C.eps*C.VP0))).^(-1);		
tempdew_fahr = 32.0+1.8.*(tempdew_kelvin-273.15d0);						
u = exp(0.1133-log(C.lambda+1d0)+0.0393.*tempdew_fahr);
t_w = 1-0.077.*(u.*m).^0.3;

%% Transmissivity after aerosol absorption
t_a = C.k_aer.^m; 

%% Transmissivity after cloud absorption / scattering
% t_cl = 1.0-0.128.*clim.C-0.346.*clim.C.^2;
t_cl = 1.0-0.233.*clim.C-0.415.*clim.C.^2;
% t_cl = 1.0-0.306.*clim.C-0.681.*clim.C.^2;

%% Incoming solar radiation
SWin = insol.TOAshade .* t_rg .* t_w .* t_a .* t_cl;

end

