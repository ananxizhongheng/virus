## python3 script
import glob
import re
import os
os.chdir('/media/wxl/Run1/wmd/gut/04host-viral/04fuzznuc_crisper_tRNA/')
# define: each GROUP must have 3 files
def getGROUP():
    G1 = [x.replace('trnaim','') for x in glob.glob("**trnaim") ]  
    G2 = [x.replace('crispr.txt','') for x in glob.glob("**crispr.txt") ] 
    G3 = [x.replace('crisprinv','') for x in glob.glob("**crisprinv") ] 
    GG = list(set(G1) & set(G2) & set(G3))
    GG.sort()
    return(GG)

GROUP = getGROUP()
print(GROUP)
f1_final = []
f2_final = []

for group in GROUP:
    # three files of each GROUP 
    file_trna = group + 'trnaim'
    dict_crispr = group + 'crispr.txt'
    file_crispr = group + 'crisprinv'
    
    # trna    
    f1 = open(file_trna,'r').readlines()
    f1_o = []
    
    for line in f1:
        if '# Sequence:' in line:
            f1_host = re.split(r'[:\s]\s*',line)[2]
        elif '#' not in line and 'NODE' in line:
            f1_virus1 = re.split(r'[+|-|:|\s]\s*',line)[4]
            f1_virus2 = re.split(r'[+|-|:|\s]\s*',line)[5]
            f1_virus = f1_virus1 + "|" + f1_virus2
            f1_out = group + " " + f1_host + " " + f1_virus
            f1_o.append(f1_out)
    f1_o = list(set(f1_o))
    f1_o.sort()
    #print(f1_o)
    
    #print("\n")
    f1_final = f1_final + f1_o
    
    
    # crispr dict
    #   this should be a seq:list(HOST_NODE) dict
    d1 = open(dict_crispr,'r').readlines()
    d1_dict = {}
    
    for line in d1:
        if 'SEQ' in line:
            d1_host =  re.split(r'[:\s]\s*',line)[1]
#        if 'Repeat' not in line and \
#                    'A' in line and \
#                    'T' in line and \
#                    'G' in line and \
#                    'C' in line:
# some sequence don't have all A,T,G,C which couldn't be located, just use '[ ]' behind 
        if '[' in line:
            d1_seq1 = re.split(r'[" "\s]\s*',line)[1]
            d1_seq2 = re.split(r'[" "\s]\s*',line)[2]
            if d1_seq1 in d1_dict:
                d1_dict[d1_seq1].append(d1_host)
            else:
                d1_dict.setdefault(d1_seq1,[]).append(d1_host)
            if d1_seq2 in d1_dict:
                d1_dict[d1_seq2].append(d1_host)
            else:
                d1_dict.setdefault(d1_seq2,[]).append(d1_host)
    #print(d1_dict["CTGTTCCCGGGCTCGAAATCCCGGGCCTCATTGAAGC"])       
    #print(d1_dict["CTCTGAAAAAAATGGAGAATGGAGAATTTG"])  
    
    #print("\n")
    
    # crispr
    f2 = open(file_crispr,'r').readlines()
    f2_o = []
    
    for line in f2:
        if '# Sequence:' in line:
            f2_virus = re.split(r'[:\s]\s*',line)[2]        
        elif '#' not in line and 'pattern' in line:
            f2_seq = re.split(r'[+|-|:|\s]\s*',line)[5]
            for f2_host in d1_dict[f2_seq]:
                f2_out = group + " " + f2_host + " " + f2_virus
                f2_o.append(f2_out)
    f2_o = list(set(f2_o))
    f2_o.sort()
    #print(f2_o)
    
    #print("\n\n")
    f2_final = f2_final + f2_o
    
f_final = f1_final + f2_final

# remove duplicates
f_final = list(set(f_final))

# sort
f_final.sort()

# write into csv file
# trna
f1_final = ["GROUP HOST VIRUS"] + f1_final
with open('./wjw-spades.tRNA.txt', 'w') as csvfile:
    for row in f1_final:
        csvfile.write(row + "\n")

# crispr
f2_final = ["GROUP HOST VIRUS"] + f2_final
with open('./wjw-spades.CRISPR.txt', 'w') as csvfile:
    for row in f2_final:
        csvfile.write(row + "\n")

# merge
f_final = ["GROUP HOST VIRUS"] + f_final
with open('./wjw-spades.merge.txt', 'w') as csvfile:
    for row in f_final:
        csvfile.write(row + "\n")



