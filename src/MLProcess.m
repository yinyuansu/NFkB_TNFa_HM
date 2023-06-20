function MLProcess(negRatio, geneSet, oriIndex, originalPatientData, outputFile)
% from patients' indices(oriIndex) of geneSet(HM/C6) and Xoriginal(_std)/yori(originalTumorcellData)
% to process ML Analysis
        % 1. [reservation] kMeansTest to get subtype (main group patients), only thinks about positive indices
        % 2. [reservation] PCA to reduce dimension (pcaRatio = 0.99; pcaOn= 1)
        % 3. svmtrain (70%)
        % 4. svmpredict (30%) (testRatio = 0.3)
% save the logs and preliminary results to outputFile

%summary: total 50 HMs(4345 genes); -TNFA HMs(4304 genes) => overlapped (157 genes in TNFA hms)

%{
geneSet = 'HM';
oriIndex = hallmark2gene('cancerGeneList.mat', 'TNFA_SIGNALING_VIA_NFKB');
%oriIndex = C6togene('cancerGeneList.mat', 'KRAS.BREAST_UP.V1_UP');
originalPatientData = '/data/blcaData_processed.mat';
outputFile = '/blca_1_log/TNFA_SIGNALING_VIA_NFKB_500.txt';
negRatio = 20;
%}


%% Initialization

% Load data to process SVM X/y array: 
load(originalPatientData); % Xoriginal(_std), yoriginal(hm), yC6(og)
if (strcmp(geneSet, 'HM'))
    yCategory = yoriginal;
elseif (strcmp(geneSet, 'C6'))
    yCategory = yC6;
end


% Get index of hallmarkgene (positive data)
pos_num = length(oriIndex); 
neg_num = negRatio*pos_num; % pos:neg

% Establish X/y
posindex = oriIndex;
rng('shuffle') %suffle the random seed
negindex = randsample(find(yCategory==-1), neg_num); %from non- hallmark/oncogenic gene to choose

% Using k-Means to get main group if necessary [reservation]
kMeansOn = 0;
if (kMeansOn)
    % Find the main group (Cancer patients) indices by K-menas
    kMgroup_num = 5;
    % Transpose matrix becuase of grouping patients
    [GroupIndex_pos, ix, B]= kMeansTest(Xoriginal_std(posindex,:)', kMgroup_num);
    mainGroupIndex_pos = find(GroupIndex_pos == ix(kMgroup_num));
else
    % Give the specific patients or total patients
    %load('ucecChoosenPatient.mat');
    mainGroupIndex_pos = (1:size(Xoriginal_std,2))';
end

y = [yCategory(posindex); yCategory(negindex)];
%X = [Xoriginal_std(posindex, mainGroupIndex_pos); Xoriginal_std(negindex, mainGroupIndex_pos)];
X = [Xoriginal_std(posindex, mainGroupIndex_pos); Xoriginal_std(negindex, mainGroupIndex_pos)];


% Using PCA to reduce dimension if necessary [reservation]
pcaRatio = 0.99;
pcaOn= 0;
[m, n] = size(X);
%  Run PCA if pcaOn==1
if (pcaOn)
    [U, S] = pca_self(X);
    k = 1; EV = diag(S);
    totalsum = sum(EV);
    while ((sum(EV(1:k))/totalsum) < pcaRatio)
        k = k+1;
    end
    Xreduce = X*U(:,1:k);
end

% Partition data for training(70%) and testing(30%)
testRatio = 0.3;
% Assign varaiables and values
posNum_test = floor(pos_num*testRatio);
negNum_test = floor(neg_num*testRatio);
if (pcaOn ==1)
    Xtest  = zeros(posNum_test+negNum_test,k);
    ytest  = zeros(posNum_test+negNum_test,1);
    Xtrain = zeros(m-size(Xtest,1),k);
    ytrain = zeros(m-size(Xtest,1),1);
else
    Xtest  = zeros(posNum_test+negNum_test,n);
    ytest  = zeros(posNum_test+negNum_test,1);
    Xtrain = zeros(m-size(Xtest,1),n);
    ytrain = zeros(m-size(Xtest,1),1);
end


% Combine pos/neg data for training/testing
pi_test  = randsample(find(y==1), posNum_test); %tumorGA_index = posindex(pi_test)
pi_train = setdiff(find(y==1), pi_test); %tumorGA_index = posindex(pi_train)
ni_test  = randsample(find(y==-1), negNum_test); %tumorGA_index = negindex(ni_test-length(posindex))
ni_train = setdiff(find(y==-1),ni_test); %tumorGA_index = negindex(ni_train-length(posindex))

if (pcaOn ==1)
    Xtest  = [Xreduce(pi_test,:); Xreduce(ni_test,:)];
    Xtrain = [Xreduce(pi_train,:); Xreduce(ni_train,:)];
else
    Xtest = [X(pi_test,:); X(ni_test,:)];
    Xtrain = [X(pi_train,:); X(ni_train,:)];
end    
ytest = [y(pi_test); y(ni_test)];
ytrain = [y(pi_train); y(ni_train)];


% Write relevant data to assigned file
fid = fopen(outputFile, 'a+'); 
if (kMeansOn)
    for i = 1:kMgroup_num
        fprintf(fid, 'Gr%d: %d(%d) (%0.2f%%):\t', i, ix(i), B(i), ...
            (B(i)/size(Xoriginal_std,2))*100);
        fprintf(fid, '%d;', find(GroupIndex_pos == ix(i))); fprintf(fid, '\n');
    end
end

fprintf(fid,'pos_train_geneIdx: ');
fprintf(fid,'%d;', posindex(pi_train)'); fprintf(fid, '\n');
fprintf(fid,'pos_test_geneIdx: ');
fprintf(fid,'%d;', posindex(pi_test)'); fprintf(fid, '\n');
fprintf(fid,'neg_train_geneIdx: ');
fprintf(fid,'%d;', negindex(ni_train-length(posindex))'); fprintf(fid, '\n');
fprintf(fid,'neg_test_geneIdx: ');
fprintf(fid,'%d;', negindex(ni_test-length(posindex))'); fprintf(fid, '\n');

if(pcaOn)
    fprintf(fid, 'kMeansOn=%d; pcaOn=%d; Before/After=%d(%0.2f%%)/%d(%0.2f%%)\n', kMeansOn, pcaOn, ...
        length(mainGroupIndex_pos), (length(mainGroupIndex_pos)/size(Xoriginal_std,2))*100, k, (k/length(mainGroupIndex_pos))*100);
end


%% ===== SVM train ===== %%
%tstart = datetime('now');
fprintf(fid, '======= SVM training start =======\n');
model = svmtrain(ytrain, sparse(Xtrain), '-s 0 -t 0 -c 1 -w1 1 -w-1 1 -m 1000');
%tend = datetime('now');
%fprintf(fid, 'The train time is: %s\n', tend-tstart);

% for linear:(w = model.SVs'*model.sv_coef, b=-model.rho) => dec_values = wX+b, 
% for nonlinear:(b=-model.rho) => dec_values = sum(model.sv_coef*k(sv, z))+b, z: test instance
fprintf(fid, 'prediction model:w= ');
w = model.SVs'*model.sv_coef;
fprintf(fid, '%f;',w); fprintf(fid, '\n');
b = -model.rho;
fprintf(fid, 'prediction model:b= %f;', b); fprintf(fid, '\n');


%% ===== SVM test ===== %%
fprintf(fid, '======= SVM predict start =======\n');
[predict_label, accuracy, dec_values] = svmpredict(ytest, sparse(Xtest), model);

fprintf(fid, 'Accuracy = %0.4f%%\n', accuracy(1));
if (~isempty(predict_label))
    tp = sum((predict_label == 1) & (ytest == 1));
    tn = sum((predict_label == -1) & (ytest==-1));
    fp = sum((predict_label == 1) & (ytest == -1));
    fn = sum((predict_label == -1) & (ytest == 1));
    prec = tp/(tp+fp); rec = tp/(tp+fn);
    F1 = (2*prec*rec)/(prec+rec);
    fprintf(fid, 'SVM_test: tp=%d; fp=%d; fn=%d\n', tp,fp,fn);
    fprintf(fid, 'SVM_test: prec = %0.4f%%; rec= %0.4f%%; F1 = %0.4f%%\n', prec*100, rec*100, F1*100);
end

% ===== all gene expression validation [pos(target gene set);neg(unexisted in any gene set)] ===== %
ni = find(yCategory==-1);
yhallmark = [yCategory(oriIndex); yCategory(ni)];
Xhallmark = [Xoriginal_std(oriIndex, :); Xoriginal_std(ni, :)];
[predict_label, accuracy, dec_values] = svmpredict(yhallmark, sparse(Xhallmark), model);

    tp = sum((predict_label == 1) & (yhallmark == 1));
    tn = sum((predict_label == -1) & (yhallmark==-1));
    fp = sum((predict_label == 1) & (yhallmark == -1));
    fn = sum((predict_label == -1) & (yhallmark == 1));
    prec = tp/(tp+fp);
    rec = tp/(tp+fn);
    F1 = (2*prec*rec)/(prec+rec);
    fprintf(fid, 'All_gene_validation: tp=%d; fp=%d; fn=%d\n', tp,fp,fn);
    fprintf(fid, 'All_gene_validation: prec = %0.4f%%; rec= %0.4f%%; F1 = %0.4f%%\n', prec*100, rec*100, F1*100);
    %[X_ROC,Y_ROC,T,AUC,OPTROCPT] = perfcurve(yhallmark, dec_values, 1);

aa = find(dec_values(length(oriIndex)+1:end)>0);
load('../matdata/cancerGeneList.mat');
cc = tumorGA(ni(aa));
load('../matdata/Gene2Oncogenic.mat');
ins = intersect(cc, Gene2ONCO(:,1));
ins_rate = length(ins)/fp;

mi = C6togene('../matdata/cancerGeneList.mat', 'MYC_UP.V1_DN', 'MYC_UP.V1_UP');
ei = C6togene('../matdata/cancerGeneList.mat', 'E2F3_UP.V1_DN', 'E2F3_UP.V1_UP');
ins_mi = intersect(cc, tumorGA(mi));
ins_ei = intersect(cc, tumorGA(ei));
fprintf(fid,'FP_geneIdx: ');
fprintf(fid,'%d;', ni(aa)'); fprintf(fid, '\n');
fprintf(fid,'FP_in_C6 = %d (%0.2f%%); FP_in_MYC = %d; FP_in_E2F3 = %d\n', length(ins), ins_rate*100, length(ins_mi), length(ins_ei));

fclose(fid);
