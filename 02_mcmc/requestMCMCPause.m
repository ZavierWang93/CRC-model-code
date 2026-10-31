function pause_file = requestMCMCPause(pathway)

pause_file = [char(pathway) '.pause'];
[fid,msg] = fopen(pause_file,'w');
if fid < 0
    error('requestMCMCPause:CreateFailed','%s',msg);
end
fprintf(fid,'Pause requested at %s\n', ...
    char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
fclose(fid);
fprintf('Pause request written to %s\n',pause_file);
end
