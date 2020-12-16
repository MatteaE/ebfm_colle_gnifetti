%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTION TO PLOT SUBSURFACE TEMPERATURE & DENSITY EVOLUTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;

%% Specify output directory & coordinates & depth of subsurface plots

outdir = '/home/ocirne94/Desktop/unifr/Tesi/code/EBFM_code_Horst/output/';
utm_x_plot = 412827;                    % specify UTM x coordinate
utm_y_plot = 5086836;                   % specify UTM y coordinate
depth_plot = 10;                        % specify depth of plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load run info
load([outdir 'runinfo.mat']);

X = grid.x_mask;                        % UTM Easting (1D)
Y = grid.y_mask;                        % UTM Northing (1D)
LX = grid.Lx; LY = grid.Ly;             % horizontal grid dimensions
L = length(X);                          % total number of glacier grid points 
dt = dtout;                             % time-step between saved data (in days)
nl = grid.nl;                           % number of vertical layers 

time_start_run = time.ts;               % time start run
time_end_run = time.te;                 % time end run

T = round((datenum(time_end_run)- ...
    datenum(time_start_run))/dt);       % total number of time-steps

%% Find nearest grid point in model output
for n=1:L
    dist = sqrt((X-utm_x_plot).^2 + (Y-utm_y_plot).^2);
    ind = find(dist==min(dist));
end

%% Read output from binary files 
T = floor((datenum(time_end_run) - datenum(time_start_run))/dt);
t = round(datenum(time_start_run)):dt:round(datenum(time_end_run));
t = t(1:T);
x = repmat(t,nl,1);

fid = fopen([outdir 'OUT_surfH.bin'],'rb');
Arr_surfH = nan(T,1);
for t=1:T
    fseek(fid,(t-1)*4*L,'bof');
    Atemp = fread(fid,L,'real*4','l');
    Arr_surfH(t,1) = Atemp(ind);
end
fclose(fid);

fid = fopen([outdir 'OUT_subdepth.bin'],'rb');
Arr_depth = nan(T,nl);
for t=1:T
    fseek(fid,(t-1)*4*L*nl,'bof');
    for n=1:nl
        Atemp = fread(fid,L,'real*4','l');
        Arr_depth(t,n) = Atemp(ind);
    end
end
fclose(fid);

surfH = (repmat(Arr_surfH,1,nl))';
y = Arr_depth' - surfH;

fid = fopen([outdir 'OUT_subD.bin'],'rb');
Arr_D = nan(T,nl);
for t=1:T
    fseek(fid,(t-1)*4*L*nl,'bof');
    for n=1:nl
        Atemp = fread(fid,L,'real*4','l');
        Arr_D(t,n) = Atemp(ind);
    end
end
fclose(fid);

fid = fopen([outdir 'OUT_subT.bin'],'rb');
Arr_T = nan(T,nl);
for t=1:T
    fseek(fid,(t-1)*4*L*nl,'bof');
    for n=1:nl
        Atemp = fread(fid,L,'real*4','l');
        Arr_T(t,n) = Atemp(ind);
    end
end
fclose(fid);


%% Figure: Density and temperature evolution
f = figure('Position',[0 0 1500 900]);

% Density subplot
ax = subplot(2,1,1);
C = (Arr_D)';
x = [x(1,:);x;x(end,:)];
y = [y(1,:)-0.00001;y;y(end,:)];
C = [C(1,:);C;C(end,:)];
C(1,:) = NaN;
C(end,:) = NaN;
xx = x(:);
yy = y(:);
zz = C(:);
F = scatteredInterpolant(xx,yy,zz,'linear');
[XX,YY] = meshgrid(min(x(:)):dt:max(x(:)),min(y(:)):0.05:max(y(:)));
G = F(XX,YY);
h = imagesc(XX(1,:),YY(:,1),G);
set(h,'alphadata',~isnan(G));
ylim([min(ylim) min(ylim)+depth_plot]);
datetickzoom('x','dd-mmm-yy'); xlabel('Date'); ylabel('Depth (m)');
c = colorbar; set(get(c,'ylabel'),'String','Density (kg m^{-3})','FontSize',12);
set(gca,'YDir','reverse');xlim([datenum(time_start_run) datenum(time_end_run)]);
set(gca,'FontSize',12);
caxis([300 1000]);

% Temperature subplot
ax = subplot(2,1,2);
C = (Arr_T)';
C = [C(1,:);C;C(end,:)];
C(1,:) = NaN;
C(end,:) = NaN;
xx = x(:);
yy = y(:);
zz = C(:)-273.15;
F = scatteredInterpolant(xx,yy,zz,'linear');
[XX,YY] = meshgrid(min(x(:)):dt:max(x(:)),min(y(:)):0.05:max(y(:)));
G = F(XX,YY);
h = imagesc(XX(1,:),YY(:,1),G); hold on;
set(h,'alphadata',~isnan(G));
ylim([min(ylim) min(ylim)+depth_plot]);
datetickzoom('x','dd-mmm-yy'); xlabel('Date'); ylabel('Depth (m)');
c = colorbar; set(get(c,'ylabel'),'String','Temperature (^{o}C)','FontSize',12);
set(gca,'YDir','reverse');xlim([datenum(time_start_run) datenum(time_end_run)]);
set(gca,'FontSize',12);

