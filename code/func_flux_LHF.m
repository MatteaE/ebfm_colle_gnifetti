%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute latent heat flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [LHF] = func_flux_LHF(C,Tsurf,clim,cond)


% 1 is glacier wind (Oerlemans and Grisogono 2002).
% 2 is standard bulk approach (Oerlemans 2000).
% 3 is bulk approach according to Suter (2002), i.e. Oke (1987).
% 4 is bulk approach according to Essery and Etchevers (2004).
LHF_scheme = 4;


if LHF_scheme == 1

    %% Turbulent exchange coefficient
    C_kat = max(C.k_turb .* (clim.T(cond)-Tsurf) .* sqrt(C.g./(C.T0 .* clim.Theta_lapse .* C.Pr)),0);
    C_turb = 0.5*(C.turb + C_kat);
    
    %% Latent heat flux (bulk equations)
    LHF = C.eps.*clim.Dair(cond).*C.Ls.*C_turb.*(clim.VP(cond)-clim.VPsat_surf)./clim.Pres(cond) .* (Tsurf<273.15) ...
    + C.eps.*clim.Dair(cond).*C.Lv.*C_turb.*(clim.VP(cond)-clim.VPsat_surf)./clim.Pres(cond) .* (Tsurf>=273.15);

elseif LHF_scheme == 2
    
    LHF = C.eps.*clim.Dair(cond).*C.Ls.*C.Ch.*clim.Wind_speed.*(clim.VP(cond)-clim.VPsat_surf(cond))./clim.Pres(cond) .* (Tsurf<273.15) ...
    + C.eps.*clim.Dair(cond).*C.Lv.*C.Ch.*clim.Wind_speed.*(clim.VP(cond)-clim.VPsat_surf(cond))./clim.Pres(cond) .* (Tsurf>=273.15);

elseif LHF_scheme == 3
    
    % Latent heat as in Suter (2002): we use the latent heat of sublimation below -10 °C surface temperature, the latent heat of vaporization at 0 °C, and a linear interpolation of the two in between.
    Lheat = C.Ls .* max(0, min(1, 0.1*(273.15-Tsurf))) ...
            + C.Lv .* max(0, min(1, 0.1*(Tsurf-263.15)));
        
    Rib_suter = C.g .* (((clim.T(cond) - Tsurf) / C.dz) / (clim.Wind_speed / C.dz)^2) ./ clim.T(cond);
    stab_fun_suter2002 = ((1 - 5 * Rib_suter).^2) .* (Rib_suter > 0) + ...
                         ((1 - 16 * Rib_suter).^0.75) .* (Rib_suter <= 0);
    
    LHF = clim.Dair(cond).*Lheat(cond).*C.k^2.*C.zu.*C.zq.*stab_fun_suter2002.* (clim.q(cond) - clim.qsat_surf(cond)) .* clim.Wind_speed ./ (C.dz ^ 2);

elseif LHF_scheme == 4
        
    % Latent heat constant formulated as in Suter (2002): we use the latent heat of sublimation below -10 °C surface temperature, the latent heat of vaporization at 0 °C, and a linear interpolation of the two in between.
    Lheat = C.Ls .* max(0, min(1, 0.1*(273.15-Tsurf))) ...
            + C.Lv .* max(0, min(1, 0.1*(Tsurf-263.15)));
    
    Rib_essery = ((C.g .* C.dz) / ((clim.Wind_speed)^2)) .* (((clim.T(cond) - Tsurf) ./ clim.T(cond)) + ((clim.q(cond) - clim.qsat_surf(cond)) ./ (clim.q(cond) + ((C.eps) / (1-C.eps)))));
    stab_fun_essery2004 = ((1 + 10 * Rib_essery).^(-1)) .* (Rib_essery > 0) + ...
                            (1 - 10 * Rib_essery .* ((1+ 10*C.Chn * sqrt(-Rib_essery)/C.fz).^(-1) )) .* (Rib_essery <= 0);
    LHF = Lheat .* (C.Chn .* stab_fun_essery2004 .* clim.Dair(cond) .* clim.Wind_speed) .* (clim.q(cond) - clim.qsat_surf(cond));
    
end

end
