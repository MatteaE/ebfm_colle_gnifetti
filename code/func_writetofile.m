function [OUT,io] = func_writetofile(OUT,io,A,grid,t,time,C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save model output to files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify variables to be written
if t==1
    io.varsout = [];
    io.unitsout = [];
    io.descout = [];
    io.dimout = [];
    if io.out_surface
        io.varsout{end+1} = 'smb';          io.unitsout{end+1} = 'm w.e. ts^{-1}';  io.descout{end+1} = 'Mass balance';
        io.varsout{end+1} = 'Tsurf';        io.unitsout{end+1} = 'K';               io.descout{end+1} = 'Surface temperature';
        io.varsout{end+1} = 'climT';        io.unitsout{end+1} = 'K';               io.descout{end+1} = 'Air temperature';
        io.varsout{end+1} = 'climP';        io.unitsout{end+1} = 'm w.e. ts^{-1}';  io.descout{end+1} = 'Precipitation';
        io.varsout{end+1} = 'climC';        io.unitsout{end+1} = 'fraction';        io.descout{end+1} = 'Cloud cover';
        io.varsout{end+1} = 'climRH';       io.unitsout{end+1} = 'fraction';        io.descout{end+1} = 'Relative humidity';
        io.varsout{end+1} = 'climPres';     io.unitsout{end+1} = 'Pa';              io.descout{end+1} = 'Air pressure';
        io.varsout{end+1} = 'climrain';     io.unitsout{end+1} = 'm w.e. ts^{-1}';  io.descout{end+1} = 'Rainfall';
        io.varsout{end+1} = 'climsnow';     io.unitsout{end+1} = 'm w.e. ts^{-1}';  io.descout{end+1} = 'Snowfall';
        io.varsout{end+1} = 'climWindspeed'; io.unitsout{end+1} = 'm s^{-1}';       io.descout{end+1} = 'Wind speed';
        io.varsout{end+1} = 'snowmass';     io.unitsout{end+1} = 'm w.e.';          io.descout{end+1} = 'Snow mass';
        io.varsout{end+1} = 'mbal';         io.unitsout{end+1} = 'm w.e.';          io.descout{end+1} = 'Cumulative mass balance';
        io.varsout{end+1} = 'mbal_stake';   io.unitsout{end+1} = 'm w.e.';          io.descout{end+1} = 'Cumulative stake mass balance';        
        io.varsout{end+1} = 'melt';         io.unitsout{end+1} = 'm w.e. ts^{-1}';  io.descout{end+1} = 'Melt';
        io.varsout{end+1} = 'refr';         io.unitsout{end+1} = 'm w.e. ts^{-1}';  io.descout{end+1} = 'Refreezing';
        io.varsout{end+1} = 'refr_seasnl_snow';    io.unitsout{end+1} = 'm w.e. ts^{-1}'; io.descout{end+1} = 'Refreezing seasonal snow';
        io.varsout{end+1} = 'refr_intacc';  io.unitsout{end+1} = 'm w.e. ts^{-1}';  io.descout{end+1} = 'Internal accumulation';
        io.varsout{end+1} = 'refr_P';       io.unitsout{end+1} = 'm w.e. ts^{-1}'; io.descout{end+1} = 'Refreezing percolating water';
        io.varsout{end+1} = 'refr_S';       io.unitsout{end+1} = 'm w.e. ts^{-1}'; io.descout{end+1} = 'Refreezing slush water';
        io.varsout{end+1} = 'refr_I';       io.unitsout{end+1} = 'm w.e. ts^{-1}'; io.descout{end+1} = 'Refreezing irreducible water';
        io.varsout{end+1} = 'runoff';       io.unitsout{end+1} = 'm w.e. ts^{-1}';  io.descout{end+1} = 'Runoff';
        io.varsout{end+1} = 'runoff_surf';  io.unitsout{end+1} = 'm w.e. ts^{-1}';  io.descout{end+1} = 'Surface runoff';
        io.varsout{end+1} = 'runoff_slush'; io.unitsout{end+1} = 'm w.e. ts^{-1}';  io.descout{end+1} = 'Slush runoff';
        io.varsout{end+1} = 'SWin';         io.unitsout{end+1} = 'W m^{-2}';        io.descout{end+1} = 'Incoming SW radiation';
        io.varsout{end+1} = 'SWout';        io.unitsout{end+1} = 'W m^{-2}';        io.descout{end+1} = 'Reflected SW radiation';
        io.varsout{end+1} = 'LWin';         io.unitsout{end+1} = 'W m^{-2}';        io.descout{end+1} = 'Incoming LW radiation';
        io.varsout{end+1} = 'LWout';        io.unitsout{end+1} = 'W m^{-2}';        io.descout{end+1} = 'Outgoing LW radiation';
        io.varsout{end+1} = 'SHF';          io.unitsout{end+1} = 'W m^{-2}';        io.descout{end+1} = 'Sensible heat flux';
        io.varsout{end+1} = 'LHF';          io.unitsout{end+1} = 'W m^{-2}';        io.descout{end+1} = 'Latent heat flux';
        io.varsout{end+1} = 'GHF';          io.unitsout{end+1} = 'W m^{-2}';        io.descout{end+1} = 'Subsurface heat flux';
        io.varsout{end+1} = 'surfH';        io.unitsout{end+1} = 'm';               io.descout{end+1} = 'Surface height';
        io.varsout{end+1} = 'alb';          io.unitsout{end+1} = 'fraction';        io.descout{end+1} = 'Albedo';   
        io.varsout{end+1} = 'climRi';       io.unitsout{end+1} = 'number';      io.descout{end+1} = 'Bulk Richardson number';
    end
    L1 = length(io.varsout);
    io.dimout(1:L1) = 2; 
    if io.out_subsurface
        io.varsout{end+1} = 'subD';         io.unitsout{end+1} = 'kg m^{-3}';       io.descout{end+1} = 'Density';
        io.varsout{end+1} = 'subT';         io.unitsout{end+1} = 'K';               io.descout{end+1} = 'Temperature';         
        io.varsout{end+1} = 'subS';         io.unitsout{end+1} = 'mm w.e.';         io.descout{end+1} = 'Slush water content'; 
        io.varsout{end+1} = 'subW';         io.unitsout{end+1} = 'mm w.e.';         io.descout{end+1} = 'Irreducible water';   
        io.varsout{end+1} = 'subdepth';     io.unitsout{end+1} = 'm';               io.descout{end+1} = 'Subsurface depth';
        io.varsout{end+1} = 'subZ';     io.unitsout{end+1} = 'm';               io.descout{end+1} = 'Layer thickness';
    end
    L2 = length(io.varsout);
    io.dimout(L1+1:L2) = 3;

%% Write grid (x, y, z, slope) to binary file
f_x = fopen([io.outdir '/OUT_x.bin'] , 'w');
fwrite(f_x, grid.x_mask,'real*4');
f_y = fopen([io.outdir '/OUT_y.bin'] , 'w');
fwrite(f_y, grid.y_mask,'real*4');
f_z = fopen([io.outdir '/OUT_z.bin'] , 'w');
fwrite(f_z, grid.z_mask,'real*4');
f_sl = fopen([io.outdir '/OUT_slope.bin'] , 'w');
fwrite(f_sl, grid.slope,'real*4');
fclose(f_x);
fclose(f_y);
fclose(f_z);
fclose(f_sl);

end


%% Update OUT.TEMP with variables to be stored
fn = io.varsout;
for i=1:numel(fn)
    temp_long = eval(['A.' fn{i}]);
    if io.dimout(i)==2
        if t>1
            OUT.TEMP.(fn{i}) = OUT.TEMP.(fn{i}) + temp_long;
        else
            OUT.TEMP.(fn{i}) = temp_long;
        end
    end
    if io.dimout(i)==3
        OUT.TEMP.(fn{i}) = temp_long;
    end
end


%% Save output to binary files
if io.save_to_binary_files
    if t==1
        for i=1:length(fn)
            io.fid(i) = fopen([io.outdir '/OUT_' fn{i} '.bin'], 'w');
        end
    end
    

    if mod(t,io.freqout)==0
        for i=1:length(fn)
            if io.dimout(i)==2
                OUT.(fn{i}) = OUT.TEMP.(fn{i})/io.freqout;
                fwrite(io.fid(i),OUT.(fn{i}),'real*4');
            end
            if io.dimout(i)==3
                OUT.(fn{i}) = OUT.TEMP.(fn{i});
                fwrite(io.fid(i),OUT.(fn{i}),'real*4');
            end
        end
    end
    
    if t==time.tn
        for i=1:length(fn)
            fclose(io.fid(i));
        end
    end
end

%% Save output to netcdf files
if io.save_to_netcdf_files
    if t==1
        delete([io.outdir '/*.nc']);
        for i=1:length(fn)
            if io.dimout(i)==2
                nccreate([io.outdir '/OUT_' fn{i} '.nc'],io.varsout{i}, ...
                    'Dimensions',{'UTM_E',grid.Ly,'UTM_N',grid.Lx,'time',Inf}, ...
                    'Format','netcdf4_classic', ...
                    'DeflateLevel',9);
            end
            if io.dimout(i)==3
                nccreate([io.outdir '/OUT_' fn{i} '.nc'],io.varsout{i}, ...
                    'Dimensions',{'UTM_E',grid.Ly,'UTM_N',grid.Lx,'layer',grid.nl,'time',Inf}, ...
                    'Format','netcdf4_classic', ...
                    'DeflateLevel',9);
                nccreate([io.outdir '/OUT_' fn{i} '.nc'],'layer','Dimensions',{'layer',grid.nl});
                ncwrite([io.outdir '/OUT_' fn{i} '.nc'],'layer',squeeze(1:grid.nl));
                ncwriteatt([io.outdir '/OUT_' fn{i} '.nc'],'layer','Units',' ');
                ncwriteatt([io.outdir '/OUT_' fn{i} '.nc'],'layer','Description','Subsurface layer number');
            end
            ncwriteatt([io.outdir '/OUT_' fn{i} '.nc'],io.varsout{i},'Units',io.unitsout{i});
            ncwriteatt([io.outdir '/OUT_' fn{i} '.nc'],io.varsout{i},'Description',io.descout{i});
            
            nccreate([io.outdir '/OUT_' fn{i} '.nc'],'UTM_N','Dimensions',{'UTM_N',grid.Lx});
            ncwrite([io.outdir '/OUT_' fn{i} '.nc'],'UTM_N',squeeze(grid.y(:,1)));
            ncwriteatt([io.outdir '/OUT_' fn{i} '.nc'],'UTM_N','Units','m');
            ncwriteatt([io.outdir '/OUT_' fn{i} '.nc'],'UTM_N','Description','UTM Northing');
            
            nccreate([io.outdir '/OUT_' fn{i} '.nc'],'UTM_E','Dimensions',{'UTM_E',grid.Ly});
            ncwrite([io.outdir '/OUT_' fn{i} '.nc'],'UTM_E',squeeze(grid.x(1,:)));
            ncwriteatt([io.outdir '/OUT_' fn{i} '.nc'],'UTM_E','Units','m');
            ncwriteatt([io.outdir '/OUT_' fn{i} '.nc'],'UTM_E','Description','UTM Easting ');
            
            nccreate([io.outdir '/OUT_' fn{i} '.nc'],'elevation','Dimensions',{'UTM_E',grid.Ly,'UTM_N',grid.Lx});
            ncwrite([io.outdir '/OUT_' fn{i} '.nc'],'elevation',grid.z');
            ncwriteatt([io.outdir '/OUT_' fn{i} '.nc'],'elevation','Units','m a.s.l.');
            ncwriteatt([io.outdir '/OUT_' fn{i} '.nc'],'elevation','Description','Surface elevation');
            
            nccreate([io.outdir '/OUT_' fn{i} '.nc'],'mask','Dimensions',{'UTM_E',grid.Ly,'UTM_N',grid.Lx});
            ncwrite([io.outdir '/OUT_' fn{i} '.nc'],'mask',grid.mask_2D');
            ncwriteatt([io.outdir '/OUT_' fn{i} '.nc'],'mask','Units','glacier (1) or no-glacier (0)');
            
            nccreate([io.outdir '/OUT_' fn{i} '.nc'],'time','Dimensions',{'time',Inf},'Format','netcdf4_classic');
            ncwriteatt([io.outdir '/OUT_' fn{i} '.nc'],'time','Units','Days since 0-Jan-0000 0:00');
        end
    end
    
    if mod(t,io.freqout)==0
        for i=1:length(fn)
            if io.dimout(i)==2
                dataOUT = nan(grid.Lx,grid.Ly);
                dataOUT(ind2sub([grid.Lx, grid.Ly],grid.ind(:))) = OUT.TEMP.(fn{i})/io.freqout;
                ncwrite([io.outdir '/OUT_' fn{i} '.nc'],io.varsout{i}, ...
                    dataOUT',[1 1 round(t/io.freqout)]);
            end
            if io.dimout(i)==3
                dataOUT = nan(grid.Lx,grid.Ly,grid.nl);
                dataOUT(ind2sub([grid.Lx, grid.Ly, grid.nl],grid.ind3(:))) = OUT.TEMP.(fn{i})/io.freqout;
                ncwrite([io.outdir '/OUT_' fn{i} '.nc'],io.varsout{i}, ...
                    permute(dataOUT,[2 1 3]),[1 1 1 round(t/io.freqout)]);
            end
            ncwrite([io.outdir '/OUT_' fn{i} '.nc'],'time',time.TCUR,round(t/io.freqout));
        end
    end
end

if mod(t,io.freqout)==0
    for i=1:length(fn)
        OUT.TEMP.(fn{i}) = 0.0;
    end
end


%% Save run info to file at end of run
if t==time.tn
    runinfo.grid = grid;
    runinfo.time = time;
    runinfo.dtout = time.dt*io.freqout;
    runinfo.varsout = io.varsout;
    runinfo.descout = io.descout;
    runinfo.unitsout = io.unitsout;
    cd(io.outdir);
    save(io.infofile,'-struct','runinfo');
    cd(io.homedir);
end

end







