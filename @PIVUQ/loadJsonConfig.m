function [cfg_data] = loadJsonConfig(configfile)
fid = fopen(configfile);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
cfg_data = jsondecode(str);
end