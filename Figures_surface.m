%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SCRIPT TO PLOT DISTRIBUTED FIELDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;

%% Specify output directory and variable to plot
%   - Open runinfo.mat to see a list of stored variables (varsout)

outdir = '/home/ocirne94/Desktop/unifr/Tesi/code/EBFM_code_Horst/output';
var = 'climP';


%% Load run info
load([outdir 'runinfo.mat']);

X = grid.x_mask;                        % UTM Easting (1D)
Y = grid.y_mask;                        % UTM Northing (1D)
UTM_x = grid.x;                         % UTM Easting (2D)
UTM_y = grid.y;                         % UTM Northing (2D) 
Lxy = [grid.Lx,grid.Ly];                % horizontal grid dimensions
L = length(X);                          % total number of glacier grid points 
dt = dtout;                             % time-step between saved data (in days)

time_start_run = time.ts;               % time start run
time_end_run = time.te;                 % time end run

ind = find(strcmp(varsout,var));        % index of selected variable ('var') in varsout

T = round((datenum(time_end_run)- ...
    datenum(time_start_run))/dt);       % total number of time-steps

%%%%%%%%%%%%%%%
%%% FIGURES
%%%%%%%%%%%%%%%

%% Example: Time-averaged spatial distribution

fid = fopen([outdir 'OUT_' var '.bin']);
Amean = zeros(L,1);
for t=1:T
    fseek(fid,(t-1)*4*L,'bof');
    Atemp = fread(fid,L,'real*4','l');
    Amean = Amean + Atemp/T;
end
fclose(fid);
temp2D = nan(Lxy(1),Lxy(2));
temp2D(ind2sub(Lxy,grid.ind(:))) = Amean(:);

figure;
pcolor(UTM_x,UTM_y,temp2D);
shading flat; xlabel('UTM x (m)'); ylabel('UTM y (m)');
axis equal; axis tight;
c = colorbar;
ylabel(c,[descout{ind} ' (' unitsout{ind} ')']);
title(descout{ind});

%% Example: Time-series spatial mean

fid = fopen([outdir 'OUT_' var '.bin']);
for t=1:T
    fseek(fid,(t-1)*4*L,'bof');
    Atemp = fread(fid,L,'real*4','l');
    temptime(t,1) = nanmean(Atemp(:));
end
fclose(fid);
tvec = datenum(time_start_run):dt:datenum(time_end_run);
tvec = tvec(1:T);

figure;
plot(tvec,temptime);
datetickzoom('x','dd-mmm-yyyy');
xlabel('Date');
ylabel([descout{ind} ' (' unitsout{ind} ')']);
title(descout{ind});