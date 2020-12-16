%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update the mass balance and snow mass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A] = func_massbalance(A,clim,C)

%% Update the mass balance
% Climatic mass balance (including internal accum. below summer surface)
A.smb = clim.snow + clim.rain - A.runoff ...
    + A.moist_deposition + A.moist_condensation ...
    - A.moist_sublimation - A.moist_evaporation;
A.mbal = A.mbal + A.smb;

% Stake mass balance (excluding internal accumulation below summer surface)
A.smb_stake = clim.snow - A.melt + A.refr_seasnl_snow ...
    + A.moist_deposition - A.moist_sublimation;
A.mbal_stake = A.mbal_stake + A.smb_stake; 

%% Update the snow mass
A.snowmass = A.snowmass + clim.snow - A.melt ...
    + A.moist_deposition - A.moist_sublimation;
A.snowmass = max(A.snowmass,0);
A.max_snowmass = sum(A.subZ.*A.subD./C.Dwater,2);
A.snowmass = min(A.snowmass,A.max_snowmass);

end


