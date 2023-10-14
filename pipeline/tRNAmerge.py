import numpy as np
import pandas as pd
import glob
import re
"""
项目需求：
整理tRNA文件到想要的形式
tRNA的代码：
aragorn -t -w -o ../06HVCs/02tRNAfromvirome/tRNA/B51-24dvftRNAw.txt B51-24.contigsdvf0.9-0.05-1.fa
aragorn -t -fons -o ../06HVCs/02tRNAfromvirome/tRNA/BHvs2tRNAfons.txt BHfinal-viral-combined.fa
注意：如果文件中出现？？？需要手调，不过这种情况极少～
"""

def getGROUP():
    G1 = [x.replace('vtRNAw.txt','') for x in glob.glob("**vtRNAw.txt") ]  
    G2 = [x.replace('vtRNAfons.txt','') for x in glob.glob("**vtRNAfons.txt") ] 
    GG = list(set(G1) & set(G2))
    GG.sort()
    return(GG)

GROUP = getGROUP()

# GROUP = ["FY","BH","B51-24","B51-11","B51","609","59","401"]
#GROUP = ["CM11-1"]
for group in GROUP:
    file_in = group + 'vtRNAw.txt'
    # data_np=np.loadtxt(file_in,dtype='str',delimiter=' ')
    data_pd = pd.read_csv(file_in, error_bad_lines=False, sep="\n",header = None)
    data_pd.columns = ['tRNA']
    data_pd.replace('\s+', '', regex=True, inplace=True)
    print(data_pd)
    data_nrows = data_pd.shape[0] # 行数
    data_ncols = data_pd.shape[1] # 列数
    print("tRNAw是 %d 行 %d 列的表格" % (data_nrows,data_ncols))

    #提取新的表格
    print("建立新表格中...")
    merge_data = pd.DataFrame()
    for i in range(data_nrows):
        if ">" in data_pd.loc[i, 'tRNA']:
            temp = i
        elif "tRNA" in data_pd.loc[i, 'tRNA']:
            # merge_data.loc[i, 1] = data_pd.loc[temp, 'tRNA']
            if "c[" in data_pd.loc[i, 'tRNA']:
                str_1 = str(data_pd.loc[i, 'tRNA']).split(r'tRNA')[1].split(r'c[')[0]
                str_2 = str(data_pd.loc[i, 'tRNA']).split(r'(')[1].split(r')')[0]
                str_3 = str(data_pd.loc[i, 'tRNA']).split(r'c[')[1].split(r']')[0]
                merge_data.loc[i, 0] =  data_pd.loc[temp, 'tRNA'] + "tRNA" + str_1+ "("+str_2 +  ")" + "c["+ str_3 + "]"
            else:
                str_1 = str(data_pd.loc[i, 'tRNA']).split(r'tRNA')[1].split(r'[')[0]
                str_2 = str(data_pd.loc[i, 'tRNA']).split(r'(')[1].split(r')')[0]
                str_3 = str(data_pd.loc[i, 'tRNA']).split(r'[')[1].split(r']')[0]
                merge_data.loc[i, 0] = data_pd.loc[
                                           temp, 'tRNA'] + "tRNA" + str_1 + "(" + str_2 + ")" + "[" + str_3 + "]"
            # print(merge_data.loc[i, 0] )

    merge_data.columns = ['tRNA+node']
    merge_data_nrows = merge_data.shape[0]  # 行数
    merge_data_ncols = merge_data.shape[1]  # 列数
    merge_data.reset_index(inplace=True, drop=True)
    print("merge_data是 %d 行 %d 列的表格" % (merge_data_nrows, merge_data_ncols))
    print(merge_data)

    file_fons = group + 'vtRNAfons.txt'
    fons_np = np.loadtxt(file_fons, dtype='str', delimiter='\t')
    print(fons_np)
    fons_data_pd = pd.DataFrame(fons_np)
    print(fons_data_pd)
    fons_data_nrows = fons_data_pd.shape[0]  # 行数
    fons_data_ncols = fons_data_pd.shape[1]  # 列数
    print("tRNAfons是 %d 行 %d 列的表格" % (fons_data_nrows, fons_data_ncols))
    fons_data_pd.columns = ['tRNA']
    merge_data_new = pd.DataFrame()

    for j in range(fons_data_nrows):
        if ">" in fons_data_pd.loc[j, 'tRNA']:
            tRNA = str(fons_data_pd.loc[j, 'tRNA']).split("tRNA-", 1)[1]
            for k in range(merge_data_nrows):
                print(merge_data.loc[k, 'tRNA+node'])
                if tRNA in merge_data.loc[k, 'tRNA+node']:
                    print(tRNA)
                    merge_data_new.loc[j, 0] = merge_data.loc[k, 'tRNA+node']
                    break
                else:
                    continue
        else:
            merge_data_new.loc[j, 0] = fons_data_pd.loc[j, 'tRNA']

    print(merge_data_new)

    file_out = group + 'vtRNAclean.txt'
    merge_data_new.to_csv(file_out, index=False, header=None)
# merge_data.to_csv("59tRNAclean.csv", index = False)


