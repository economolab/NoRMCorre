function [chans, planes] = parseTiffChansPlanesInfo(info)
id = info.ImageDescription;

tok = regexp(id, 'channels=(\d*)', 'tokens');
if isempty(tok)
    disp('Could not read number of channels in tif.  Assuming 1 channel');
    chans = 1;
else
    chans = str2double(tok{1}{1});
end

tok = regexp(id, 'slices=(\d*)', 'tokens');

if isempty(tok)
    tok = regexp(id, 'frames=(\d*)', 'tokens');
end
if isempty(tok)
    disp('Could not read number of z planes in tif.  Assuming 1 plane');
    planes = 1;
else
    planes = str2double(tok{1}{1});
end