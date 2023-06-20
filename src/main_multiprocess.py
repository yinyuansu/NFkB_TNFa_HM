import os, sys, time
from multiprocessing import Pool

def external_call(i):
	os.system('matlab /minimize /nosplash /nodesktop /r \"main '+i)

cancers = ["blca", "brca", "chol", "coad", "esca", "hnsc", "kich", "kirc", "kirp", "lihc", "luad", "lusc", "prad", "stad", "thca", "ucec"]
ex_queue =  []
path = 'D:/NFkB_TNFa_HM/'
subfolder = 'log_cHM_120_1000/'
negRatio = '20'
for cancer in cancers:
	for i in range(1,3):
		#give parameters "path+subfolder+logfolder+dataName+geneSet" to main.m
		item = path+' '+subfolder+' '+cancer+'_'+str(i)+'_log '+cancer+' HM '+negRatio+'\"'
		#print(item)
		ex_queue.append(item)
#print(ex_queue)

if __name__ == '__main__':
	with Pool(16) as p:
		print('16 processes were started at ' + str(time.asctime()) + ' !')
		p.map(external_call, ex_queue)
