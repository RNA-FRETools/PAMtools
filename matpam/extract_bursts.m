% Extract photon bursts binned in PAM and exported to a Matlab figure

function extract_bursts

fig = gcf;
axObjs = fig.Children;

if strcmp(fig.Name, 'PAM: PIE Analysis with Matlab')
    disp('No figure to extract.')
    disp('You need to export the burst trace from PAM first')
    return
end

ax = axObjs(2).Children;
data.time = ax.XData;
data.counts = ax.YData;


% save file
[file, path] = uiputfile('*.json', 'Save JSON file', 'bursts.json');
if file
    fileID = fopen([path file],'wt');
    jsontext = jsonencode(data);
    fprintf(fileID,'%s',jsontext);
    fclose(fileID);
end
close(fig)
end
