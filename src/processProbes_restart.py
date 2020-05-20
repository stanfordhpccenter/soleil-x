import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# SET INPUTS
size = '512by256'
date = '20200319'
folder_1 = date + '_' + size + '/'
folders_2 = ["1/", "2/", "3/", "9/", "10/"]
folder_3 = "job1/"
folder_4 = "sample0/"
pathfile_out = 'results_' + date + '/' 


column_names = "0", "1", "2", "3", "4", "5"
summary = pd.DataFrame(columns= column_names)

for sample in folders_2:
    foldername = folder_1 + sample + folder_3 + folder_4
    for i in range(6):
        filename = 'probe' + str(i) + '.csv'
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
summary

os.makedirs(pathfile_out)
summary.to_csv(pathfile_out + size + '_summary.csv')
