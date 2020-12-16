%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ENERGY BALANCE - FIRN MODEL (EBFM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% References: 
%%%   * Van Pelt et al. (2012), doi:10.5194/tc-6-641-2012
%%%   * Van Pelt and Kohler (2015), doi:10.3189/2015JoG14J223
%%%   * Van Pelt et al. (2019), doi:10.5194/tc-2019-53
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;

%% Model setup
[grid,time,io]              = func_init_params();
[C]                         = func_init_constants();
[grid]                      = func_init_grid(grid,io);
[A,clim,insol,OUT]          = func_init_arrays(C,grid,io);

% To log the exchange coefficient.
% A.logvar = zeros(time.tn,1);
% A.logvar2 = zeros(time.tn,1);


%% Time-loop
for t=1:time.tn
    
    %% Print time to screen
    [time] = func_printtime(t,time);

    %% Read and prepare climate input
    [clim,A] = func_loadclimate(C,grid,clim,io,t,time,A);
    
    %% Surface energy balance model
    [A,insol] = func_energybalance(C,A,clim,insol,t,time,grid);
    
    %% Snow/firn model
    [A] = func_snowmodel(C,A,clim,time.dt,grid,time);
    
    %% Mass balance
    [A] = func_massbalance(A,clim,C);
    
    %% Write output to files
    [OUT,io] = func_writetofile(OUT,io,A,grid,t,time,C);
        
end

%% Save restart-files
func_createbootfile(A,io);
