% Extract 2D histogram drawn in PAM and exported to a Matlab figure

% Specify any additional dictionary elements to be included in the json
% output file as name-value pairs. 
% 
% Currently implemented names are:
% - photons_per_window (int)
% - crosstalk (float)
% - direct_excitation (float)
% - gamma_factor (float)
%
% Example: extract_2Dplot('photons_per_window', 5, 'direct_excitation', 0.046, ... 
%                         'crosstalk', 0.11, 'gamma_factor', 0.89)
%

function extract_2Dplot(varargin)

fig = gcf;
axObjs = fig.Children;
try
    dataObjs = axObjs.Children;
catch 
    disp('No figure to extract.')
    disp('You need to export the 2D plot from BurstExplorer first')
    return
end

for i = 1:length(dataObjs)
    ax = dataObjs(i).Children;
    if strcmp(dataObjs(i).Tag, 'Axes_1D_X')
        barAxis = ax(strcmp(get(ax, 'type'), 'bar'));
        data.X1D.x = barAxis.XData;
        data.X1D.y = barAxis.YData;
        lineAxis = ax(strcmp(get(ax, 'type'), 'line'));
        data.X1DfitComp.x = {};
        data.X1DfitComp.y = {};
        data.limits.x = dataObjs(i).XLim;
        k = 1;
        for j=1:6
            if strcmp(lineAxis(j).LineStyle, '-')
                data.X1DfitSum.x = lineAxis(j).XData;
                data.X1DfitSum.y = lineAxis(j).YData;
            else
                if length(lineAxis(j).XData) > 1
                    data.X1DfitComp.x{k} = lineAxis(j).XData;
                    data.X1DfitComp.y{k} = lineAxis(j).YData;
                    k = k+1;
                end
            end
        end
    elseif strcmp(dataObjs(i).Tag, 'Axes_1D_Y')
        barAxis = ax(strcmp(get(ax, 'type'), 'bar'));
        data.Y1D.x = barAxis.XData;
        data.Y1D.y = barAxis.YData;
        lineAxis = ax(strcmp(get(ax, 'type'), 'line'));
        data.Y1DfitComp.x = {};
        data.Y1DfitComp.y = {};
        data.limits.y = dataObjs(i).XLim;
        k = 1;
        for j=1:6
            if strcmp(lineAxis(j).LineStyle, '-')
                data.Y1DfitSum.x = lineAxis(j).XData;
                data.Y1DfitSum.y = lineAxis(j).YData;
            else
                if length(lineAxis(j).XData) > 1
                    data.Y1DfitComp.x{k} = lineAxis(j).XData;
                    data.Y1DfitComp.y{k} = lineAxis(j).YData;
                    k = k+1;
                end
            end
        end
    elseif strcmp(dataObjs(i).Tag, 'Colorbar')
        data.cmapLimits = dataObjs(i).Limits;
    elseif strcmp(dataObjs(i).Tag, 'Axes_General')
        
        isContour = strcmp(get(ax, 'type'), 'contour');
        if any(isContour)
            contourAxis = ax(isContour);
            contourAxis = contourAxis(end);
            data.contour.x = contourAxis.XData;
            data.contour.y = contourAxis.YData;
            data.contour.z = contourAxis.ZData;
            data.contour.levels = contourAxis.LevelList;
        else
            data.contour.x = string(nan);
            data.contour.y = string(nan);
            data.contour.z = string(nan);
            data.contour.levels = string(nan);
        end
        
        isImage = strcmp(get(ax, 'type'), 'image');
        if any(isImage)
            imageAxis = ax(isImage);
            imageAxis = imageAxis(end);
            data.image.x = imageAxis.XData;
            data.image.y = imageAxis.YData;
            data.image.z = imageAxis.CData;
        else
            data.image.x = string(nan);
            data.image.y = string(nan);
            data.image.z = string(nan);
        end
        
        isScatter = strcmp(get(ax, 'type'), 'scatter');
        if any(isScatter)
            scatterAxis = ax(isScatter);
            scatterAxis = scatterAxis(end);
            data.scatter.x = scatterAxis.XData;
            data.scatter.y = scatterAxis.YData;
            if length(data.scatter.x) <= 2
                scatterNull = 1;
            else
                scatterNull = 0;
            end
        else
            scatterNull = 1;
        end
        if scatterNull
            data.scatter.x = string(nan);
            data.scatter.y = string(nan);
        end
        
        isPatch = strcmp(get(ax, 'type'), 'patch');
        if any(isPatch)
            patchAxis = ax(isPatch);
            patchAxis = patchAxis(end);
            data.hex.x = patchAxis.XData;
            data.hex.y = patchAxis.YData;
            data.hex.z = patchAxis.CData;
        else
            data.hex.x = string(nan);
            data.hex.y = string(nan);
            data.hex.z = string(nan);
        end
         
    end
end

p = inputParser;
default_photons_per_window = nan;
addParameter(p,'photons_per_window',default_photons_per_window,@(x) assert(isnumeric(x) && isscalar(x) && mod(x,1) == 0 && x>0, 'Value must be a positive integer'));
addParameter(p,'gamma_factor',default_photons_per_window,@(x) assert(isnumeric(x) && isscalar(x) && x>0, 'Value must be a positive float'));
addParameter(p,'crosstalk',default_photons_per_window,@(x) assert(isnumeric(x) && isscalar(x) && x>0, 'Value must be a positive float'));
addParameter(p,'direct_excitation',default_photons_per_window,@(x) assert(isnumeric(x) && isscalar(x) && x>0, 'Value must be a positive float'));

parse(p,varargin{:});
f = fieldnames(p.Results);
for i = 1:length(f)
    if ~isnan(p.Results.(f{i}))
        data.(f{i}) = p.Results.(f{i});
    end
end 

% save file
[file, path] = uiputfile('*.json', 'Save JSON file', 'data.json');
if file
    fileID = fopen([path file],'wt');
    jsontext = jsonencode(data);
    fprintf(fileID,'%s',jsontext);
    fclose(fileID);
end
close(fig)
end

