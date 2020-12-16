%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute subsurface heat flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [GHF] = func_flux_GHF(Tsurf,A,GHF_C,cond)

%% Thermal heat flux between surface and subsurface
GHF = GHF_C(cond).*(A.subT(cond,3)-Tsurf(:));

end