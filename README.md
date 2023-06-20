## Ensemble learning model for identifying the hallmark genes of NFκB/TNF signaling pathway in cancers
**We have developed a powerful computational framework that facilitates the identification of genes associated with specific pathways and carcinogenesis. In this demonstration, we applied our framework to analyze gene expression profiles obtained from patients across 16 distinct cancer types. By utilizing the ensemble learning model, it is able to identify genes that contribute to both disease development and the NFκB/TNF signaling pathway. This approach allows to uncover the crucial genes that play a role in the pathogenesis of various cancers. Furthermore, our framework has the capability of conducting pan-cancer inference by integrating multiple ensemble models derived from diverse cancer types. This integration ensures a comprehensive understanding of the underlying mechanisms across different cancer types.**

##

### Quick start guide
We provide comprehensive lists that include 50 hallmark gene sets, 189 oncogenic signature gene sets, and gene expression profiles of 16 cancer types for further analysis and application.
Here's a quick start guide to get you started:	
* 1. Specify the gene sets: Open the "Hallmark_list.txt" file and specify the gene sets you are interested in analyzing.
* 2. Gene expression profiles: Make sure you have the standardized gene expression profiles of patients for the cancer types. These profiles should be available in the "data" folder.
* 3. Run the ensemble learning models: Execute the "main_multiprocess.py" script to parallelly train ensemble learning models for all the cancer types. This step will generate the necessary log files for further analysis.
* 4. Parse coefficients: Run the "getSVMpar_byHallmarks.m" script to parse the coefficients of all SVM models. This step will extract the relevant information from the trained models.
* 5. Extract vote counts: Finally, execute the "extractPredictedPositive.m" script to extract the vote counts of the whole genome for each cancer type. This step will provide you with valuable insights into the predicted positive genes associated with the specific cancer types.


### Code and data

`src/:`
* Hallmark_list.txt: A text file containing a list of target hallmark gene sets.
* Hallmark_list_all.txt: A text file containing a list of all 50 hallmark gene sets from MSigDB.
* oncogenic_category.txt: A text file containing a list of all 189 oncogenic gene sets from MSigDB.
* hallmark2gene.m: A script that finds all the genes for given hallmark gene sets.
* C6togene.m: A script that finds all the genes for given oncogenic gene sets.
* main_multiprocess.py: A script that handles the multiprocess for the main.m process.
* main.m: The main process that handles the MLProcess.m implementation.
* MLProcess.m: A script that implements the ensemble learning models and saves logs to the log_cHM_120_1000 folder.
* getSVMpar_byHallmarks.m: A script that parses all the coefficients of SVM models and saves the structure to SVMpar.mat.
* extractPredictedPositive.m: A script that extracts predicted positive genes from ensemble learning models and saves the results to wgpredict_bymodels.mat.
* kMeansTest.m: A script for performing k-means clustering.
* pca_self.m: A script for PCA dimension reduction.

`matdata/:`
* cancerGeneList.mat: A file containing lists of whole genome and non-hallmark genes in patients.
* Gene2CancerHallmarks.mat: A file containing a list of whole hallmark genes (geneID, geneset#, list of hallmark gene sets).
* Gene2Oncogenic.mat: A file containing a list of whole oncogenic signature genes (geneID, geneset#, list of oncogenic gene sets).
* SVMpar.mat: A file containing the coefficients of all SVM models in the ensemble learning, extracted from getSVMpar_byHallmarks.m.
* wgpredict_bymodels.mat: A file containing vote counts of whole genome by cancers from ensemble learning models.
* wgpredict_bymodels.txt: A text file format of wgpredict_bymodels.mat for further use.

`data/:` Gene expression profiles and gene set labels for all cancer types.

`log_cHM_120_1000/:` Log files of ensemble learning models.

`Third-party software:` The SVM algorithm of our ensemble learning model is implemented using the libSVM library. We have provided the executable binary files for convenience. For more details, please refer to the *libSVM website* and the *src/libSVM_readme* file.
    
    src/:
	    svmpredict.mexw64: Executable binary file for the SVM prediction function.
	    svmtrain.mexw64: Executable binary file for the SVM training function.
	    libSVM_readme: Readme file with additional information about the libSVM library.

## Citation

