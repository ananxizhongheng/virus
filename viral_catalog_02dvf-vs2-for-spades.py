#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 10:19:15 2022
项目要求：spades组装的contig获得*vs2+dvf0.9-0.05merge.csv结果

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
#print(GROUP)

#GROUP = ["F3F4"]
for group in GROUP:
    # 读VS2
    file_in = group + 'viralseqname.txt'
    file_out = group + 'uniquename.csv'
    vs2_pd = pd.read_csv(file_in, error_bad_lines=False, sep="\t")
    vs2_pd_nrows = vs2_pd.shape[0]
    # print(vs2_pd)


    # 读dvf0.7/0.9
    # file7_in = "/Users/liuxinwu/Desktop/dvf/"+ group + ".contigs" + "0.7-0.05-1kdvf.txt"
    file9_in = "../04posdvf/" + group + ".contigs" + "0.9-0.05-10kdvf.txt"
    # dvf7_pd = pd.read_csv(file7_in, error_bad_lines=False, sep=",")
    # dvf7_pd_nrows = dvf7_pd.shape[0]
    dvf9_pd = pd.read_csv(file9_in, error_bad_lines=False, sep="\t")
    dvf9_pd_nrows = dvf9_pd.shape[0]
    # dvf7_pd.columns = ['contig name','dvf_length','dvf_score>0.7','dvf_pvalue<0.05','dvf_qvalue']
    dvf9_pd.columns = ['contig name', 'dvf_length', 'dvf_score>0.9', 'dvf_pvalue<0.05']

    # 统计大于length的值
    length = 10000
    dvf9_count_len = 0
    for i in range(dvf9_pd_nrows):
        if dvf9_pd.loc[i, 'dvf_length'] > length:
            dvf9_count_len = dvf9_count_len + 1
    # print(dvf9_pd)
    print("%s 的DVF0.9包含 %d 个contig (其中length大于 %d 的有 %d 个)" % (group, dvf9_pd_nrows,length,dvf9_count_len ))

    # 整理一下VS2的表格
    merge_data = pd.DataFrame()

    for i in range(vs2_pd_nrows):
        merge_data.loc[i,'contig name'] = str(vs2_pd.loc[i,'seqname']).split(r'|')[0]
        merge_data.loc[i, 'vs2_length'] = str(vs2_pd.loc[i,'seqname']).split(r'|')[0].split(r'length_')[1].split(r'_cov_')[0]
        merge_data.loc[i, 'vs2'] = str(vs2_pd.loc[i,'seqname']).split(r'|')[2]

    # print(merge_data)
    # merge_data.set_index('contig name', drop=False)
    merge_data = merge_data.drop_duplicates(['contig name'], keep='first')
    merge_data.to_csv(file_out, index=False)
    merge_data = pd.read_csv(file_out, error_bad_lines=False, sep=",")
    vs2_pd_uniquenrows = merge_data.shape[0]
    # print(vs2_pd_uniquenrows)
    # print(merge_data)
    vs2_count_len = 0
    for i in range(vs2_pd_uniquenrows):
        # print(merge_data.loc[i, 'vs2_length'])
        if float(merge_data.loc[i, 'vs2_length']) > length:
            vs2_count_len = vs2_count_len + 1
    print("%s 的VS2包含 %d 个contig (其中length大于 %d 的有 %d 个)" % (group, vs2_pd_uniquenrows, length, vs2_count_len))



    # 求交集/并集
    # merge7_data = pd.merge(merge_data, dvf7_pd, how='outer', on='contig name')
    # both7_data = pd.merge(merge_data, dvf7_pd, how='inner', on='contig name')
    merge9_data = pd.merge(merge_data, dvf9_pd, how='outer', on='contig name')
    # both9_data = pd.merge(merge_data, dvf9_pd, how='inner', on='contig name')



    # merge7_data_nrows = merge7_data.shape[0]
    # both7_data_nrows = both7_data.shape[0]
    merge9_data_nrows = merge9_data.shape[0]
    # both9_data_nrows = both9_data.shape[0]
    # for i in range(merge7_data_nrows):
    #     if pd.isnull(merge7_data.at[i, 'vs2_length']):
    #         merge7_data.loc[i, 'vs2_length'] = merge7_data.loc[i, 'dvf_length']
    # merge7_data = merge7_data.rename(columns={'vs2_length': 'length'})
    # merge7_data = merge7_data.drop(['dvf_length'], axis=1)

    for i in range(merge9_data_nrows):
        if pd.isnull(merge9_data.at[i, 'vs2_length']):
            merge9_data.loc[i, 'vs2_length'] = merge9_data.loc[i, 'dvf_length']
    merge9_data = merge9_data.rename(columns={'vs2_length': 'length'})
    merge9_data = merge9_data.drop(['dvf_length'], axis=1)

    # 统计大于length的值
    # merge7_count_len = 0
    # for i in range(merge7_data_nrows):
    #     if int(merge7_data.loc[i, 'length']) > length:
    #         merge7_count_len = merge7_count_len + 1
    merge9_count_len = 0
    for i in range(merge9_data_nrows):
        if int(merge9_data.loc[i, 'length']) > length:
            merge9_count_len = merge9_count_len + 1
    # both7_count_len = 0
    # for i in range(both7_data_nrows):
    #     if int(both7_data.loc[i, 'vs2_length']) > length:
    #         both7_count_len = both7_count_len + 1
    # both9_count_len = 0
    # for i in range(both9_data_nrows):
    #     if int(both9_data.loc[i, 'vs2_length']) > length:
    #         both9_count_len = both9_count_len + 1

    # 保存
    # file7_both_out= "/Users/liuxinwu/Desktop/vs2+dvf0.7/"+ group + "vs2+dvf0.7-0.05both.csv"
    # file7_merge_out= "/Users/liuxinwu/Desktop/vs2+dvf0.7/" + group + "vs2+dvf0.7-0.05merge.csv"
    # file9_both_out = "/Users/liuxinwu/Desktop/vs2+dvf0.9/" + group + "vs2+dvf0.9-0.05both.csv"
    file9_merge_out = "../07dvf_vs2/" + group + "vs2+dvf0.9-0.05merge.csv"


    # merge7_data.to_csv(file7_merge_out, index=False)
    # both7_data.to_csv(file7_both_out, index=False)
    merge9_data.to_csv(file9_merge_out, index=False)
    # both9_data.to_csv(file9_both_out, index=False)

    # print("%s 的VS2-DVF0.7merge文件包含 %d 个contig(其中length大于 %d 的有 %d 个)" % (group, merge7_data_nrows, length, merge7_count_len))
    # print("%s 的VS2-DVF0.7both文件包含 %d 个contig(其中length大于 %d 的有 %d 个)" % (group, both7_data_nrows, length, both7_count_len))
    print("%s 的VS2-DVF0.9merge文件包含 %d 个contig(其中length大于 %d 的有 %d 个)" % (group, merge9_data_nrows, length, merge9_count_len))
    # print("%s 的VS2-DVF0.9both文件包含 %d 个contig(其中length大于 %d 的有 %d 个)" % (group, both9_data_nrows, length, both9_count_len))


