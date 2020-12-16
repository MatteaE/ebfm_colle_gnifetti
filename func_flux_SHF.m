%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute sensible heat flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SHF] = func_flux_SHF(C,Tsurf,Theta_surf,clim,cond)


% 1 is glacier wind (Oerlemans and Grisogono 2002).
% 2 is standard bulk approach (Oerlemans 2000).
% 3 is bulk approach according to Suter (2002), i.e. Oke (1987).
% 4 is bulk approach according to Essery and Etchevers (2004).
SHF_scheme = 4;

if SHF_scheme == 1

    %% Turbulent exchange coefficient
    C_kat = max(C.k_turb .* (clim.T(cond)-Tsurf) .* sqrt(C.g./(C.T0 .* clim.Theta_lapse .* C.Pr)),0);
    C_turb = 0.5*(C.turb + C_kat);
    
    %% Sensible heat flux (bulk equations)
    SHF = clim.Dair(cond).*C.Cp.*C_turb.*(clim.T(cond)-Tsurf);
    
elseif SHF_scheme == 2
    
    SHF = clim.Dair(cond).*C.Cp.*C.Ch.*clim.Wind_speed.*(clim.T(cond)-Tsurf);

elseif SHF_scheme == 3
  
    Rib_suter = C.g .* (((clim.T(cond) - Tsurf) / C.dz) / (clim.Wind_speed / C.dz)^2) ./ clim.T(cond);
    stab_fun_suter2002 = ((1 - 5 * Rib_suter).^2) .* (Rib_suter > 0) + ...
                         ((1 - 16 * Rib_suter).^0.75) .* (Rib_suter <= 0);
                     
    SHF = clim.Dair(cond).*C.Cp.*C.k^2.*C.zu.*C.zt.*stab_fun_suter2002.*clim.Wind_speed.*(clim.T(cond) - Tsurf) ./ (C.dz ^ 2);

elseif SHF_scheme == 4
    
    Rib_essery = ((C.g .* C.dz) / ((clim.Wind_speed)^2)) .* (((clim.T(cond) - Tsurf) ./ clim.T(cond)) + ((clim.q(cond) - clim.qsat_surf(cond)) ./ (clim.q(cond) + ((C.eps) / (1-C.eps)))));
    stab_fun_essery2004 = ((1 + 10 * Rib_essery).^(-1)) .* (Rib_essery > 0) + ...
                            (1 - 10 * Rib_essery .* ((1+ 10*C.Chn * sqrt(-Rib_essery)/C.fz).^(-1) )) .* (Rib_essery <= 0);
    
    SHF = (C.Chn .* stab_fun_essery2004 .* clim.Dair(cond) .* C.Cp .* clim.Wind_speed) .* (clim.Theta(cond) - Theta_surf);
    
end

end
