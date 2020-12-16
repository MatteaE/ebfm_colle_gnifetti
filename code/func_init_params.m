%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [grid, time, io] = func_init_params()
clearvars;

%% Time parameters
time.ts = '1-Jan-2014 00:00';                                            % date start run
time.TS = datenum(time.ts); 
time.te = '31-Dec-2015 23:00';                                          % date end run
time.TE = datenum(time.te);
time.ti = '1-Jan-2003 00:00';                                            % date start climate record
time.TI = datenum(time.ti);
% time.dt = 0.125;                                               % timestep (days)
time.dt = 0.0416666666666666666666667;                        % hourly series                  
time.tn = round((time.TE-time.TS)/time.dt)+1;                           % nr of timesteps


%% Grid parameters
grid.utmzone = 32;                                                      % UTM zone
grid.max_subZ = 0.10;                                                    % maximum layer thickness (m)
grid.nl = 250;                                                          % number of vertical layers
grid.doubledepth = 0;                                                   % double vertical layer depth at layer nr defined in grid.split (1=yes, 0=no)
grid.split = [11,21,26];                                                % vertical layer nr at which layer depth doubles


%% Input/output parameters
io.homedir = '/home/user/EBFM_paper_final';      % home directory
io.rebootdir = [io.homedir '/reboot/'];                                 % reboot directory
io.outdir = [io.homedir '/output/'];                                    % output directory

% 20 m grids
io.gridfile = [io.homedir '/grid/grid_20m.mat'];                            % grid file containing DEM + mask
io.climfile = [io.homedir '/climate/summerbias/climate_1h.mat'];                      % climate input file
io.clim_accum_file = [io.homedir '/accumulation/clim/accumulation_climatology_20m.tif']; % annual accumulation climatology grid


io.out_surface = 1;                                                     % write surface variables to files (1=yes, 0=no)
io.out_subsurface = 1;                                                  % write subsurface variables to files (1=yes, 0=no)
io.freqout = 1;                                                         % frequency of storing output (every n-th time-step)

io.readbootfile = 0;                                                    % read initial conditions from file (1=yes, 0=no)
io.bootfilein = 'boot_file_20m_1h_mat';                                     % bootfile to be read
io.writebootfile = 1;                                                   % write file for rebooting (1=yes, 0=no)  
io.bootfileout = 'boot_file_20m_1h_mat';                                   % bootfile to be written

io.infofile = 'runinfo.mat';                                            % file to store run information

io.save_to_binary_files = 1;                                            % save output to binary files (1=yes, 0=no)
io.save_to_netcdf_files = 0;                                            % save output to netcdf files (1=yes, 0=no)

end

