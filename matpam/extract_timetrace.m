% Extract photon bursts binned in PAM and exported to a Matlab figure
%
% Specify optional start and end time limits for the exported trace as
% name-value pairs
% 
% - start_time (float)
% - end_time (float)
% 
% Example: extract_timetrace('start_time', 0.001, 'end_time', 1)

function extract_bursts(varargin)

fig = gcf;
axObjs = fig.Children;

if strcmp(fig.Name, 'PAM: PIE Analysis with Matlab')
    disp('No figure to extract.')
    disp('You need to export the intensity time trace from PAM first')
    return
end

ax = axObjs(2).Children;

p = inputParser;
default_start = ax.XData(1);
default_end = ax.XData(end);
addParameter(p,'start_time',default_start,@(x) assert(x>=default_start && x<default_end, sprintf('Value must be a within %s and %s', default_start, default_end)));
addParameter(p,'end_time',default_end,@(x) assert(x>default_start && x<=default_end, sprintf('Value must be a within %s and %s', default_start, default_end)));
parse(p,varargin{:});

[~, start_ind] = min(abs(p.Results.start_time-ax.XData));
[~, end_ind] = min(abs(p.Results.end_time-ax.XData));

data.time = ax.XData(start_ind:end_ind);
data.counts = ax.YData(start_ind:end_ind);

% save file
[file, path] = uiputfile('*.json', 'Save JSON file', 'timetrace.json');
if file
    fileID = fopen([path file],'wt');
    jsontext = jsonencode(data);
    fprintf(fileID,'%s',jsontext);
    fclose(fileID);
end
close(fig)
end
