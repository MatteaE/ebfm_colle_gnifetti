%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute reflected shortwave radiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [SWout,A] = func_flux_SWout(C,time,A,SWin)

%% Knap and Oerlemans + Bougamont albedo formulation (see Van Pelt et al. 2019)
soil_cond = A.subSOIL(:,1)==1;

A.tstar(A.Tsurf==273.15) = C.tstar_wet;
A.tstar(A.Tsurf<273.15) = C.tstar_dry + min(273.15-A.Tsurf(A.Tsurf<273.15),10)*C.tstar_K;

A.alb_snow = A.alb_snow - (A.alb_snow-C.alb_firn)./A.tstar*time.dt; 
A.alb_snow(A.timelastsnow==time.TCUR) = C.alb_fresh;
A.alb_snow(A.snowmass==0) = C.alb_fresh;

A.alb(~soil_cond,1) = A.alb_snow(~soil_cond,1) + (C.alb_ice-A.alb_snow(~soil_cond,1)).*exp(-A.snowmass(~soil_cond,1).*1d3./C.dstar);
A.alb(soil_cond,1) = C.soil_albedo;


%% Reflected SW radiation
SWout = SWin .* A.alb;



end

