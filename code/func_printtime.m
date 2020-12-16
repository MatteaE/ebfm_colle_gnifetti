%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Print simulation year and day number to screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time] = func_printtime(t,time)

%% Print time to screen while running
time.TCUR = time.TS+(t-1)*time.dt;
tempdate = datevec(time.TCUR);

fprintf('%6s','Year: ');    fprintf('%3i',tempdate(1));
fprintf('%9s','   Month: ');  fprintf('%3i',tempdate(2));
fprintf('%7s','   Day: ');  fprintf('%3i',tempdate(3));
fprintf('%8s','   Hour: ');  fprintf('%3i',tempdate(4));
fprintf('%6s','   TS: '); fprintf('%5i\n',t);

end



