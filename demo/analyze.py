import xlrd
import pandas as pd
import os
import numpy as np

print(1)
xl_workbook = xlrd.open_workbook("3_lib_final.xlsx")
sheet_names = xl_workbook.sheet_names()
xl_sheet = xl_workbook.sheet_by_name(sheet_names[0])

prom_lib = {}
ref_bc = []
ref_bc_last10 = []
ref_p = []

bc_dict = {}
prom_dict = {}

#column3 = barcode (reverse complemented to make 5'to3' sequence), column4 = promoter sequence without ATG+Barcode
for rnum in range(1,xl_sheet.nrows):
    prom_lib[xl_sheet.cell(rnum,5).value] = xl_sheet.cell(rnum,3).value
    ref_bc.append(xl_sheet.cell(rnum,5).value)
    ref_bc_last10.append(xl_sheet.cell(rnum,3).value[-10:])
    ref_p.append(xl_sheet.cell(rnum,3).value)
    ref_dict = dict(zip(ref_bc_last10, ref_bc))

def gc(s):
    return float((s.count("G") + s.count("C"))/len(s))   

def read_dict(readfile):
    dict = {}
    file = open(readfile, "r")
    for line in file.readlines():
        key = line.split(", ")[0]
        value = float(line.split(", ")[1])
        dict[key] = value
    return dict


#Parse RNA reads

masterdf = pd.DataFrame.from_dict(prom_lib,orient="index")
masterdf.columns = ["promoter"]
masterdf["GC_content"] = [gc(x) for x in masterdf["promoter"]]

head = []
Rsum = []
RNA_files = []
c = 0
for files in os.listdir("0_bccounts/"):
    
    if "RNA" in files:
        c+=1
        if c<10:
            head.append("0" + str(c))
        else:
            head.append((str(c)))
        RNA_files.append("0_bccounts/" + files)
        
for pos in range(len(head)):
    tempRNA = pd.DataFrame(index = prom_lib.keys()) 
    tempRNA[head[pos] + "_RNAcounts"] = pd.DataFrame.from_dict(read_dict(RNA_files[pos]),orient="index")
    masterdf = pd.concat([masterdf,tempRNA],axis = 1)

for pos in range(len(head)):
    RNAsum = (masterdf[head[pos] + "_RNAcounts"]).sum()

    Rsum.append(RNAsum)

    normRNAs = pd.Series(masterdf[head[pos] + "_RNAcounts"]/RNAsum)
    normRNAs.columns = [head[pos] + "_RNAnorm"]
    masterdf[head[pos] + "_RNAnorm"] = normRNAs

masterdf.to_csv("RNA_parse.csv")

#Parse DNA reads

masterdf = pd.DataFrame.from_dict(prom_lib,orient="index")
masterdf.columns = ["promoter"]
masterdf["GC_content"] = [gc(x) for x in masterdf["promoter"]]

head = []
Dsum = []
DNA_files = []
c = 0
for files in os.listdir("0_bccounts/"):
    
    if "DNA" in files:
        c+=1
        if c<10:
            head.append("0" + str(c))
        else:
            head.append((str(c)))
        DNA_files.append("0_bccounts/" + files)
        
for pos in range(len(head)):
    tempDNA = pd.DataFrame(index = prom_lib.keys()) 
    tempDNA[head[pos] + "_DNAcounts"] = pd.DataFrame.from_dict(read_dict(DNA_files[pos]),orient="index")
    masterdf = pd.concat([masterdf,tempDNA],axis = 1)

for pos in range(len(head)):
    DNAsum = (masterdf[head[pos] + "_DNAcounts"]).sum()
    Dsum.append(DNAsum)

    normDNAs = pd.Series(masterdf[head[pos] + "_DNAcounts"]/DNAsum)
    normDNAs.columns = [head[pos] + "_DNAnorm"]
    masterdf[head[pos] + "_DNAnorm"] = normDNAs

masterdf.to_csv("DNA_parse.csv")


DNAdf = pd.read_csv("DNA_parse.csv",index_col = [0])
RNAdf = pd.read_csv("RNA_parse.csv",index_col = [0])

Dcutoff = 10
Rcutoff = 10

masterdf = pd.DataFrame.from_dict(prom_lib,orient="index")
masterdf.columns = ["promoter"]
masterdf["GC_content"] = [gc(x) for x in masterdf["promoter"]]

head = []
c = 0
for files in os.listdir("0_bccounts/"):
    if "DNA" in files:
        c+=1
        if c<10:
            head.append("0"+str(c))
        else:
            head.append(str(c))
#basically, divide each DNA norm col by RNA norm col, except for the ones thats less than cutoff counts

#first set up the dataframe

for pos in range(len(head)):
    masterdf = pd.concat([masterdf,DNAdf[head[pos] + "_DNAcounts"],RNAdf[head[pos] + "_RNAcounts"]],axis = 1)

for pos in range(len(head)):
    DNAsum = (masterdf[head[pos] + "_DNAcounts"]).sum()
    RNAsum = (masterdf[head[pos] + "_RNAcounts"]).sum()

    print(DNAsum)
    print(RNAsum)
    print()
    normDNAs = pd.Series(masterdf[head[pos] + "_DNAcounts"]/DNAsum)
    normDNAs.columns = [head[pos] + "_DNAnorm"]
    normRNAs = pd.Series(masterdf[head[pos] + "_RNAcounts"]/RNAsum)
    normRNAs.columns = [head[pos] + "_RNAnorm"]
    masterdf[head[pos] + "_DNAnorm"] = normDNAs
    masterdf[head[pos] + "_RNAnorm"] = normRNAs
    act = pd.Series((np.log10(normRNAs[p]/normDNAs[p]) if masterdf[head[pos] + "_RNAcounts"][p] >= Rcutoff and masterdf[head[pos] + "_DNAcounts"][p] >= Dcutoff else np.nan for p in range(len(masterdf.index))),name = head[pos] + "_activ",index = prom_lib.keys())
    masterdf[head[pos] + "_activ"] = act 
    
masterdf.to_csv("master_data.csv")