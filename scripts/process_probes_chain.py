import pandas as pd
import numpy as np
import os
import errno

folder = '/usr/workspace/adcock4/soleilx/src/20200319_512by256/'
pathfile_out = '/p/gpfs1/adcock4/results_stochasticBC/512by256/'
samples = ['1', '2', '3', '9', '10']

try:
    os.makedirs(pathfile_out)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

column_names = "0", "1", "2", "3", "4", "5"
summary = pd.DataFrame(columns= column_names)

for sample in samples:
    foldername = folder + sample
    for i in range(6):
        filename_start = '/job0/sample0/probe' + str(i) + '.csv'
        filepath_start = foldername + filename_start
        data_start = pd.read_csv(filepath_start, delimiter='\t')
        T_initial = data_start.at[0, 'AvgFluidT']

        filename_end = '/job1/sample0/probe' + str(i) + '.csv'
        filepath_end = foldername + filename_end
        data_end = pd.read_csv(filepath_end, delimiter='\t')
        N_rows_end = data_end.shape[0] -1 #subtract off header
        T_final = data_end.at[N_rows_end, 'AvgFluidT']
        T_diff = T_final - T_initial
        summary.at[sample,str(i)]=T_diff

summary.at['max']=summary.max(axis=0)
summary.at['min']=summary.min(axis=0)
summary.at['mean']=summary.mean(axis=0)
summary.at['variance']=summary.var(axis=0)
summary.at['(max-min)/mean (%)']= \
    (summary.loc['max']-summary.loc['min'])/summary.loc['mean']*100


summary.to_csv(pathfile_out + 'summary_966418_964713_965161_965183_965186.csv')
