#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 13:08:09 2022
项目要求：spades组装的contig获得both,onlydvf，onlyvs2结果

@author: anan
"""

import pandas as pd
import os
import re
import glob

##获得目前的位置
os.chdir('/Users/anan/Documents/北京大学/油藏/油藏实验/ORIGC/结果/wmd/06vs2name/')
os.getcwd()
def getGROUP():
    GG = [x.replace('viralseqname.txt','') for x in glob.glob("**viralseqname.txt") ]
    GG.sort()
    return(GG)

GROUP = getGROUP()

#GROUP = ["BH4"]


for group in GROUP:
    # 读VS2
    file_in = group + 'vs2name.txt'
    file_out = group + 'uniquename.csv'
    vs2_pd = pd.read_csv(file_in, error_bad_lines=False, sep="\t", header = None)
    vs2_pd_nrows = vs2_pd.shape[0]

    # 整理一下VS2的表格
    vs2_pd.columns = ['contig name']
    for i in range(vs2_pd_nrows):
        name = str(vs2_pd.loc[i,'contig name'])
        vs2_pd.loc[i,'contig name'] = name.split(r'|')[0].split(r'vs2_')[1]
        vs2_pd.loc[i, 'vs2tag name'] = name.split(r'>')[1]
        vs2_pd.loc[i, 'vs2_length'] = name.split(r'|')[0].split(r'length_')[1].split(r'_cov_')[0]
        vs2_pd.loc[i, 'vs2'] = "vs2-" + name.split(r'|')[2]
    # print(vs2_pd)
    vs2_count_len = 0
    length = 10000
    for i in range(vs2_pd_nrows):
        if int(vs2_pd.loc[i, 'vs2_length']) >= length:
            vs2_count_len = vs2_count_len + 1
    print("%s 的VS2包含 %d 个contig (其中length不小于 %d 的有 %d 个)" % (group, vs2_pd_nrows, length, vs2_count_len))


    # 读dvf0.9
    file9_in = "../04posdvf/" + group + ".contigs0.9-0.05-10kdvf_name.txt"
    dvf9_pd = pd.read_csv(file9_in, error_bad_lines=False, sep=",",header=None)
    dvf9_pd_nrows = dvf9_pd.shape[0]
    dvf9_pd.columns = ['contig name']
    for i in range(dvf9_pd_nrows):
        dvf9_pd.loc[i,'contig name'] = dvf9_pd.loc[i,'contig name']
        dvf9_pd.loc[i, 'dvftag name'] = group + "dvf_" + dvf9_pd.loc[i,'contig name']
        dvf9_pd.loc[i, 'dvf'] = "dvfscore0.9"
        dvf9_pd.loc[i, 'dvf_length'] = str(dvf9_pd.loc[i,'contig name']).split(r'length_')[1].split(r'_cov_')[0]
    # print(dvf9_pd)
    # 统计大于length的值
    length = 10000
    dvf9_count_len = 0
    for i in range(dvf9_pd_nrows):
        if int(dvf9_pd.loc[i, 'dvf_length']) >= length:
            dvf9_count_len = dvf9_count_len + 1
    print("%s 的DVF0.9包含 %d 个contig (其中length不小于 %d 的有 %d 个)" % (group, dvf9_pd_nrows,length,dvf9_count_len ))





    # 求交集/并集
    merge9_data = pd.merge(vs2_pd, dvf9_pd, how='outer', on='contig name')
    both9_data = pd.merge(vs2_pd, dvf9_pd, how='inner', on='contig name')
    # print(merge9_data)
    # print(both9_data)
    merge9_data_nrows = merge9_data.shape[0]
    both9_data_nrows = both9_data.shape[0]

    # 保存
    file9_both_out = "../07dvf_vs2/" + group + "vs2+dvf0.9-0.05both.csv"
    file9_merge_out = "../07dvf_vs2/" + group + "vs2+dvf0.9-0.05merge.csv"
    for i in range(merge9_data_nrows):
        if pd.isnull(merge9_data.at[i, 'vs2_length']):
            merge9_data.loc[i, 'vs2_length'] = merge9_data.loc[i, 'dvf_length']
    merge9_data = merge9_data.rename(columns={'vs2_length': 'length'})
    merge9_data = merge9_data.drop(['dvf_length'], axis=1)
    both9_data = both9_data.rename(columns={'vs2_length': 'length'})
    both9_data = both9_data.drop(['dvf_length'], axis=1)
    merge9_data.to_csv(file9_merge_out, index=False)
    both9_data.to_csv(file9_both_out, index=False)

    # 提取txt
    vs2_data = pd.DataFrame()
    dvf_data = pd.DataFrame()
    both_data = pd.DataFrame()
    vs2_data_out = "../07dvf_vs2/" + group + "onlyvs2seqkit-10k.txt"
    dvf_data_out = "../07dvf_vs2/" + group + "onlydvfseqkit-10k.txt"
    both_data_out = "../07dvf_vs2/" + group + "bothseqkit-10k.txt"
    vs2_num = 0
    dvf_num = 0
    both_num = 0
    length = 10000
    merge9_count_len = 0
    both9_count_len = 0
    for i in range(merge9_data_nrows):
        if int(merge9_data.loc[i, 'length']) >= length:
            merge9_count_len = merge9_count_len + 1
            if pd.isnull(merge9_data.at[i, 'dvftag name']):
                vs2_data.loc[vs2_num, 'contig name'] = merge9_data.loc[i, 'vs2tag name']
                vs2_num = vs2_num + 1
            else:
                if pd.isnull(merge9_data.at[i, 'vs2tag name']):
                    dvf_data.loc[dvf_num, 'contig name'] = merge9_data.loc[i, 'contig name']
                    dvf_num = dvf_num + 1
                else:
                    both_data.loc[both_num, 'contig name'] = merge9_data.loc[i, 'contig name']
                    both_num = both_num + 1
    for j in range(both9_data_nrows):
        if int(both9_data.loc[j, 'length']) >= length:
            both9_count_len = both9_count_len + 1

    dvf_data.to_csv(dvf_data_out, index=False, header = False)
    vs2_data.to_csv(vs2_data_out, index=False, header = False)
    both_data.to_csv(both_data_out, index=False, header=False)
    print("%s 的VS2-DVF0.9merge文件包含 %d 个contig(其中length不小于 %d 的有 %d 个)" % (group, merge9_data_nrows, length, merge9_count_len))
    print("%s 的VS2-DVF0.9both文件包含 %d 个contig(其中length不小于 %d 的有 %d 个)" % (group, both9_data_nrows, length, both9_count_len))



