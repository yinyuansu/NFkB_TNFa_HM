%ver distcomp
function main(path, subfolder, logfolder, dataName, geneSet, negRatio)


close all; clc;

%{
path = 'D:/NFkB_TNFa_HM/';
subfolder = 'log_cHM_120_1000/';
logfolder = 'ucec_1_log'; %output folder
dataName = 'ucec'; %processed data
geneSet = 'HM'; %using gene set (HM,C6,...)
negRatio = 20; %negnum = negRatio*posnum
%}
mkdir(strcat(path, subfolder, logfolder));

% Get the Hallmark_gene/oncogenic_gene/... set list
if (strcmp(geneSet, 'HM'))
    CList = textread(strcat(path,'src/Hallmark_list.txt'),'%s');
elseif (strcmp(geneSet, 'C6'))
    CList = textread(strcat(path,'src/oncogenic_category.txt'),'%s');
end

patientGA = strcat(path,'matdata/cancerGeneList.mat'); %tumorGA is the same for all cancer
for i = 1:length(CList)
    % Get patients' indices of hallmark_gene/oncogenic_gene/Cancer_gene (positive data)
    if (strcmp(geneSet, 'HM'))
        oriIndex = hallmark2gene(patientGA, CList{i});
    elseif (strcmp(geneSet, 'C6'))
        oriIndex = C6togene(patientGA, CList{i});
    end
    % e.g. hallmark2gene('cancerGeneList.mat', 'ESTROGEN_RESPONSE_LATE', 'ESTROGEN_RESPONSE_EARLY');
    outputFile = strcat(path, subfolder, logfolder, '/', CList{i}, '_500.txt'); % needed assigned output file
    fid = fopen(outputFile, 'a+'); 
    
    tic
    %parfor (j = 1:1000, 6)	% only use 6 cores
    for (j = 1:500)	
        fprintf(fid, '\n:################: %s (%s_%d) :################:\n', CList{i}, dataName, j);
        MLProcess(str2double(negRatio), geneSet, oriIndex, strcat(path,'data/',dataName,'Data_processed.mat'), outputFile); % the only needed processed data
    end
    time2 = toc;
    
    fprintf(fid, 'Total run time of 500*%s = %g\n', CList{i}, time2);
    fclose(fid);
    %parpool close
    %delete(gcp('nocreate'));
   
end
