import os as os
import pandas as pds    
    
wd_path = "/home/nick/git/split_test_FEM/"
os.chdir(wd_path)
sample_list = os.listdir("./samples_var/") #NOTE THE sample data folder here
num_samples = len(sample_list)

for ii in range(0,num_samples):
    sample_list[ii] = ("./samples_var/"+sample_list[ii])

headers = ["sample_type", "theta", "cutp", "cutn", "seed_num","mes", "stress_mean", "sample_num"]

df = pds.read_csv(sample_list[0],names=headers)

for ii in range(len(df.sample_type)):
    df.sample_type[ii] = df.sample_type[ii].split("/")[-1].split(".")[0]
    df.sample_num[ii] = 0

for ii in range(1,num_samples):
    df_tmp = pds.read_csv(sample_list[ii],names=headers)
    for jj in range(len(df_tmp.sample_type)):
        df_tmp.sample_type[jj] = df_tmp.sample_type[jj].split("/")[-1].split(".")[0]
        df_tmp.sample_num[jj] = ii
    df = pds.concat([df,df_tmp], ignore_index=True)

df.to_csv("sample_data_var.csv", sep=",")  #output file name in main dir
