% Parsing "weight" and "bias" of SVM from given hallmarks and cancer types, 
% and save coefficients matrix of 1000 models to SVMpar.mat
% In each {cancer} file, (row: model#; column: 1~end-1 is weight, end is bias).


hallmarks = textread('Hallmark_list.txt','%s');
cancers = ["blca", "brca", "chol", "coad", "esca", "hnsc", "kich", "kirc", "kirp", "lihc", "luad", "lusc", "prad", "stad", "thca", "ucec"];
log_path = '../log_cHM_120_1000/';
matdata_path = '../matdata/';
run_times = 500;
file_num = 2;

% Parsing partition (w,b)
pars = ["prediction model:w= ", "prediction model:b= "];
SVMpar = struct();
% parsing parameters(w,b) of SVM to SVMpar.mat (SVMpar.{hm}.cancers)
matFile = strcat(matdata_path, 'SVMpar.mat');

for hm = 1:length(hallmarks)
    for cn = 1:length(cancers) % parsing file by each cancer
        PAR = [];
        for fn = 1:file_num
            parsedFile = strcat(log_path, cancers(cn), '_', num2str(fn), '_log\', hallmarks{hm}, '_500.txt');
            fid = fopen(parsedFile);
            while 1
                tline = fgetl(fid);
                if ~ischar(tline), break, end
                % parsing weight of each feature(w)
                if ~isempty(findstr(tline, pars(1)))
                    [token,remain] = strtok(tline, '='); % Cut the header
                    w = str2num(strrep(remain, '= ', '')); % get w
                    PAR = [PAR;w];
                end
                % parsing bias of each model(b)
                if ~isempty(findstr(tline, pars(2)))
                    [token,remain] = strtok(tline, '='); % Cut the header
                    b = str2num(strrep(remain, '= ', '')); % get b
                    PAR = [PAR;b];
                end
            end
            fclose(fid);
        end
        RAR = reshape(PAR,[],run_times*file_num);
        SVMpar.(hallmarks{hm}).(cancers(cn)) = RAR';
    end
    save(matFile, 'SVMpar');
end
