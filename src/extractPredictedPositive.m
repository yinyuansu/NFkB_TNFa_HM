% Extract predicted positive genes (hm) from 1000 training models in each cancer type (SVMpar.mat),
% and save results to 'wgpredict_bymodels.mat'

cancers = ["blca", "brca", "chol", "coad", "esca", "hnsc", "kich", "kirc", "kirp", "lihc", "luad", "lusc", "prad", "stad", "thca", "ucec"];
log_path = '../log_cHM_120_1000/';
matdata_path = '../matdata/';
hallmarks = textread('Hallmark_list.txt','%s');

wgPredictbymodels = struct();

for hm = 1:length(hallmarks)
    matFile = strcat(matdata_path, 'SVMpar.mat');
    load (matFile);

    % predict whole genome data by 1000 models and save to "wgPredictbymodels" structure
    load(strcat(matdata_path,'cancerGeneList.mat'));
    wgPredictbymodels_tmp = tumorGA;
    oriIndex = hallmark2gene(strcat(matdata_path,'cancerGeneList.mat'), hallmarks{hm}); %find hm_gene
    load('../data/blcaData_processed.mat');
    % column2: hm_status(0:ni_gene, 1:other_hm,  2:hm_gene)
    ni     = find(yoriginal == -1);
    all_hm = setdiff(1:length(yoriginal),ni); other_hm = setdiff(all_hm,oriIndex);
    wgPredictbymodels_tmp(:,2) = "1"; % 1:other_hm gene
    wgPredictbymodels_tmp(ni,2) = "0"; % 0:ni_gene
    wgPredictbymodels_tmp(oriIndex,2) = "2"; % 2:hm_gene
    % column3: C6_status(0/1)
    C6_idx = find(yC6 == 1);
    wgPredictbymodels_tmp(:,3) = "0";
    wgPredictbymodels_tmp(C6_idx,3) = "1";
    wgPredictbymodels_tmp = ["gene", "hm_status", "C6"; wgPredictbymodels_tmp];

    %% main function
    for cn = 1:length(cancers)
        % get Xoriginal_std and yoriginal
        originalPatientData = strcat('../data/', cancers(cn),'Data_processed.mat');
        load(originalPatientData);
        %oriIndex = C6togene('cancerGeneList.mat', 'MYC_UP.V1_DN');
        oriIndex = hallmark2gene(strcat(matdata_path,'cancerGeneList.mat'), hallmarks{hm});
        model_par = SVMpar.(hallmarks{hm}).(cancers(cn));
        ni = find(yoriginal == -1);
        %ni = find(yC6 == -1);

        % predict positive genes by 1000 models
        %model_predict = (model_par(:,1:end-1)*Xoriginal_std(oriIndex,:)') + model_par(:,end);
        % predict whole genome data by 1000 models
        model_predict = (Xoriginal_std*model_par(:,1:end-1)') + model_par(:,end)';
        model_vote = sum((model_predict>0),2);
        model_vote = [cancers(cn);model_vote];
        wgPredictbymodels_tmp = [wgPredictbymodels_tmp,model_vote];
    end
    wgPredictbymodels.(hallmarks{hm}) = wgPredictbymodels_tmp;
    save(strcat(matdata_path, 'wgpredict_bymodels.mat'), 'wgPredictbymodels');
    
    % copy "wgPredictbymodels_tmp" to wgPredict_bymodels({hm}).txt file
    outputFile = strcat(matdata_path,'wgPredict_bymodels(', hallmarks{hm},').txt');
    fid = fopen(outputFile, 'a+'); 
    for raw = 1:length(wgPredictbymodels_tmp)
        for column = 1:size(wgPredictbymodels_tmp,2)
            fprintf(fid,'%s\t', wgPredictbymodels_tmp(raw,column));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
end        
