import pandas as pd
import numpy as np
import os
import errno

folder = '/p/gpfs1/adcock4/'
pathfile_out = '/p/gpfs1/adcock4/results_stochasticBC/128by64/'
#samples = ['961602', '961603', '961604', '961605', '961606']
samples = ['961607', '961608', '961613', '961610', '961611']

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
        filename = '/sample0/probe' + str(i) + '.csv'
        filepath = foldername + filename
        data = pd.read_csv(filepath, delimiter='\t')
        N_rows = data.shape[0] -1 #subtract off header
        T_initial = data.at[0, 'AvgFluidT']
        T_final = data.at[N_rows, 'AvgFluidT']
        T_diff = T_final - T_initial
        summary.at[sample,str(i)]=T_diff

summary.at['max']=summary.max(axis=0)
summary.at['min']=summary.min(axis=0)
summary.at['mean']=summary.mean(axis=0)
summary.at['variance']=summary.var(axis=0)
summary.at['(max-min)/mean (%)']= \
    (summary.loc['max']-summary.loc['min'])/summary.loc['mean']*100


summary.to_csv(pathfile_out + 'summary_961607-8_10-11_13.csv')
