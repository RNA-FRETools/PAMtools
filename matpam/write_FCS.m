% FCS data from workspace to file
%
% FCS : struct
%       columns are "time", "data", "error", "fit" and "residuals"
% average_counts : column vector (n x 1 double) => manually copy from table into a column vector [c1;c2;...] on the Matlab command line


function write_FCS(FCS, average_counts)

%---------------------------------
% fitmodel definitions
FCS_D_trip_nogamma.func = '(1/P(1))*(1+((P(6)/(1-P(6))*exp(-x/P(5)/1e-6)))).*(1./(1+4*(P(2)*1e-12)*x/(P(3)*1e-6)^2)).*(1./sqrt(1+4*(P(2)*1e-12)*x/(P(4)*1e-6)^2))+P(7)';
FCS_D_trip_nogamma.parameters = {'N', 'D', 'w_xy', 'w_z', 'tau_t', 'A_t', 'y0'};
fitmodels.FCS_D_trip_nogamma = FCS_D_trip_nogamma;

FCS_D_2trip_nogamma.func = '(1/P(1))*(1+((P(6)/(1-P(6))*exp(-x/P(5)/1e-6)))).*(1+((P(8)/(1-P(8))*exp(-x/P(7)/1e-6)))).*(1./(1+4*(P(2)*1e-12)*x/(P(3)*1e-6)^2)).*(1./sqrt(1+4*(P(2)*1e-12)*x/(P(4)*1e-6)^2))';
FCS_D_2trip_nogamma.parameters = {'N', 'D', 'w_xy', 'w_z', 'tau_t', 'A_t', 'tau_t2', 'A_t2', 'y0'};
fitmodels.FCS_D_2trip_nogamma = FCS_D_2trip_nogamma;
%---------------------------------

% get fit model
model_function_input = strrep(char(FCS.Model),'@(P,x)','');
n_models = numel(fieldnames(fitmodels));
model_names = fieldnames(fitmodels);
model_found = 0;
for i = 1:n_models
    model_name_sele = model_names{i};
    model_function_sele = fitmodels.(model_names{i}).func;
    model_parameters_sele = fitmodels.(model_names{i}).parameters;
    if strcmp(model_function_sele, model_function_input)
        fprintf('-> model found!\n')
        model_found = 1;
        break
    end
end
if ~model_found
    fprintf('-> model does not exist\n');
    return
end
    

    
% rewrite model with parameters
model_var_sele = model_function_input;
for p = 1:numel(fitmodels.(model_name_sele).parameters)
    model_var_sele = strrep(model_var_sele, sprintf('P(%d)', p), char(fitmodels.(model_name_sele).parameters(p)));
end


for f = 1:size(FCS.Params,1)
    if exist('filepath', 'var')
        [file,filepath] = uiputfile(sprintf('%s_fcs.txt',fullfile(filepath,FCS.FileName{f})));
    else
        [file,filepath] = uiputfile(sprintf('%s_fcs.txt',FCS.FileName{f}));
    end
    
    % data and fits
    fileID = fopen(fullfile(filepath,file),'w');
    if fileID ~= -1
        fprintf(fileID,'%s\t%s\t%s\t%s\t%s\n',FCS.Graphs{1}{:});
        fprintf(fileID, '%0.10f\t%0.10f\t%0.10f\t%0.10f\t%0.10f\n',FCS.Graphs{f+1}');
        fclose(fileID);
        fprintf('-> data and fits written to file %s\n', file)
        
        
        % fit parameters            
        n = length(model_parameters_sele);
        if exist('average_counts', 'var')
            n = n+1;
        end
        file_param = strrep(file,'.txt','_param.txt');
        fileID = fopen(fullfile(filepath,file_param),'w');
        fprintf(fileID,'# %s\n',model_name_sele);
        fprintf(fileID,'# %s\n',model_var_sele);
        formatSpec = repmat(['%s\t'; '%f\t'], 1,n);
        formatSpec(:,end) = 'n';
        
        if exist('average_counts', 'var')
            fprintf(fileID,formatSpec(1,:), model_parameters_sele{:}, 'average_counts');
            fprintf(fileID,formatSpec(2,:), FCS.Params(f,:),average_counts(f));
        else
            fprintf(fileID,formatSpec(1,:), model_parameters_sele{:});
            fprintf(fileID,formatSpec(2,:),FCS.Params(f,:));
        end
        fprintf(fileID,'\n');
        fclose(fileID);
        fprintf('-> model parameters written to file %s\n', file_param)
    end

end




