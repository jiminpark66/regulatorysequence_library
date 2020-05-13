import xlrd
import os
import gzip
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
from datetime import datetime
from Bio import SeqIO
from collections import Counter
import pandas as pd
from collections import defaultdict
import math

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
    
DNAfilepaths = []
RNAfilepaths = []
fp = "merged"

for subdir, dirs, files in os.walk(fp):
    for file in files:
        if "fastq" in file and "gz" in file:
            a = subdir.split("/")[-1] + "/"
            b = file
            path = os.path.join(a,b) 
            if "D" in path:
                DNAfilepaths.append(path)
            elif "-" in path:
                RNAfilepaths.append(path)
            else:
                raise("errrrrrrrr")

def exp_err(phred_arr):
    E = 0.0
    for q in phred_arr:
        E += math.pow(10,(-q/10))
    return E   

def tabulate(key, counts):
    # add 1 to count for this key, or add it to the dictionary    
    if key in counts:
        counts[key] += 1
    else:
        counts[key] = 1
        
def tabulate_read(key,value,dict):
    # add value to count for this key, or add it to the dictionary    
    if key in dict:
        dict[key].append(value)
    else:
        dict[key] = [value]
        
def write_dict(dict, outputfile):
    # write sequences and counts to a text file
    print (outputfile)
    file = open(outputfile, 'w+')
    for w in sorted(dict, key=dict.get, reverse=True):
        file.write('{seq}, {num}\n'.format(seq=w, num=dict[w]))
    file.close()

# DNA parsing
def parse_dna(file):
    fadapter = "GGATCC" #BamHI
    radapter = "CTGCAG" #PstI

    #alignment scoring
    match = 0
    mismatch = -1
    gap_open = -1
    gap_extend = -1
    score_dist = []

    sref_bc = set(ref_bc)
    exp_err_dist = []
    read_trim_len_dist = []
    prom_len_dist = []
    bc_len_dist = []
    
    count = 0
    both_adaptercount = 0
    
    lowq_count = 0
    missingadapter_count = 0
    frag_count = 0
    badbc_count = 0
    goodbc_perfectalign_count = 0
    goodbc_goodalign_count = 0
    goodbc_badalign_count = 0
    
    lowq = defaultdict(int)
    missingadapter = defaultdict(int)
    frag = defaultdict(int)
    badbc = defaultdict(int)
    goodbc_perfectalign = defaultdict(list)
    goodbc_goodalign = defaultdict(list)
    goodbc_badalign = defaultdict(list)
    
    goodbc_perfectalign_bccounts = defaultdict(int)
    goodbc_goodalign_bccounts = defaultdict(int)
    
    print("Parsing of " + file + " started at:" + str(datetime.now()) + "\n")

    handle = gzip.open(file,"rt")
    for rec in SeqIO.parse(handle, 'fastq'):
        count +=1
        qscore = rec.letter_annotations["phred_quality"]
        read = str(rec.seq)
        error = exp_err(qscore)
        exp_err_dist.append(error)

        if error <3:
            
            if fadapter in read and radapter in read:
                both_adaptercount +=1
                               
                fpos = read.find(fadapter)
                rpos = read.find(radapter)
                read_trim = read[fpos+len(fadapter):rpos]
                read_trim_len_dist.append(len(read_trim))
# modified                   
                if 117 >= len(read_trim) >= 113:
                    prom = read_trim[:-15]
                    prom_len_dist.append(len(prom))
                    bc = read_trim[-12:]
                    bc_len_dist.append(len(bc))
                    
# modified                  
                    sbc = str(bc)
                    if len(bc)==12 and 102>=len(prom)>=98:
                        
                        if sbc in sref_bc:
                            if prom == prom_lib[sbc]:
                                goodbc_perfectalign_count += 1
                                tabulate(sbc,goodbc_perfectalign_bccounts)
                                tabulate_read(sbc,prom,goodbc_perfectalign)
                                score_dist.append(0)
                            else:
                                score = pairwise2.align.globalms(prom,prom_lib[sbc],match,mismatch,gap_open,gap_extend,score_only = True)
                                score_dist.append(score)
# modified
                                if score >= (math.ceil(len(prom)*0.04)*-1):
                                    goodbc_goodalign_count +=1
                                    tabulate(sbc,goodbc_goodalign_bccounts)
                                    tabulate_read(sbc,prom,goodbc_goodalign)
                                else:
                                    goodbc_badalign_count += 1
                                    tabulate_read(sbc,prom,goodbc_badalign)
                        else:
                            badbc_count +=1
                            tabulate(read_trim,badbc)
                                    
                    else:
                        badbc_count +=1
                        tabulate(read_trim,badbc)
                else:
                    frag_count +=1  
                    tabulate(read_trim,frag)
                    
            else:
                missingadapter_count += 1
                tabulate(read,missingadapter)
            
        else:
            tabulate(read,lowq)
            lowq_count += 1
            
        #if count >100000:
            #break
            
    combined = {}
    combined.update(goodbc_perfectalign_bccounts)
    combined.update(goodbc_goodalign_bccounts)
        
    combined_bccounts = {}
    for key in combined:
        if key in goodbc_goodalign_bccounts and key in goodbc_perfectalign_bccounts:
            combined_bccounts[key] = goodbc_goodalign_bccounts[key] + goodbc_perfectalign_bccounts[key] 
        elif key in goodbc_goodalign_bccounts:
            combined_bccounts[key] = goodbc_goodalign_bccounts[key]
        else:
            combined_bccounts[key] = goodbc_perfectalign_bccounts[key]
            
    print("total count: " + str(count))
    print("")
    print("lowq count: " +str(lowq_count) + " = " + "{0:.2f}".format(float(lowq_count/count)))
    print("missing_adapter count: "+ str(missingadapter_count) + " = " + "{0:.2f}".format(float(missingadapter_count/count)))
    print("frag count: " + str(frag_count) + " = " + "{0:.2f}".format(float(frag_count/count)))
    print("badbc count: " + str(badbc_count) + " = " + "{0:.2f}".format(float(badbc_count/count)))
    print("goodbc_badalignment count: " + str(goodbc_badalign_count) + " = " + "{0:.2f}".format(float(goodbc_badalign_count/count)))
    print("goodbc_goodalignment count: " + str(goodbc_goodalign_count) + " = " + "{0:.2f}".format(float(goodbc_goodalign_count/count)))
    print("goodbc_perfectalignment count: " + str(goodbc_perfectalign_count) + " = " + "{0:.2f}".format(float(goodbc_perfectalign_count/count)))
    print("")

    filename = file.split("/")[0] + "/" + file.split("/")[-1].split(".")[0]
    runname = file.split("/")[0]
    fname = file.split("/")[-1].split(".")[0]
    d = pd.DataFrame(list(combined_bccounts.items()), columns = ['Barcode', 'Counts'])
    d = d.set_index('Barcode')
    d.to_csv("05_bc_counts/" + fname + ".csv")
    d2 = pd.read_csv("05_bc_counts/" + fname + ".csv")

    BC = ""
    BC_lst = []
    for BC_rev in d2['Barcode']:
        BC = Seq(BC_rev).reverse_complement()
        BC = ''.join(BC)
        BC_lst.append(BC)
    d2['Barcode'] = BC_lst
    d2.to_csv("05_bc_counts/" + fname + ".csv", index=False)
    
# modified
    #write to files
    write_dict(combined_bccounts,"0_bccounts/" + fname + ".txt")
    write_dict(lowq,"1_lowq/" + fname + ".txt")
    write_dict(missingadapter,"2_missingadapter/" + fname + ".txt")
    write_dict(frag,"3_frag/" + fname + ".txt")
    write_dict(badbc,"4_badbc/" + fname + ".txt")
    write_dict(goodbc_badalign,"5_goodbc_badalign/" + fname + ".txt")
    write_dict(goodbc_goodalign,"6_goodbc_goodalign/" + fname + ".txt")
    write_dict(goodbc_perfectalign,"7_goodbc_perfectalign/" + fname + ".txt")

    write_dict(goodbc_goodalign_bccounts,"8_goodbc_goodalign_bccounts/" + fname + ".txt")
    write_dict(goodbc_perfectalign_bccounts,"9_goodbc_perfectalign_bccounts/" + fname + ".txt")
    
    with open("10_log_files/" + fname + ".txt","w+") as f:
        f.write("total count: " + str(count) + "\n")
        f.write("" + "\n")
        f.write("lowq count: " +str(lowq_count) + " = " + "{0:.2f}".format(float(lowq_count/count)) + "\n")
        f.write("missing_adapter count: "+ str(missingadapter_count) + " = " + "{0:.2f}".format(float(missingadapter_count/count)) + "\n")
        f.write("frag count: " + str(frag_count) + " = " + "{0:.2f}".format(float(frag_count/count)) + "\n")
        f.write("badbc count: " + str(badbc_count) + " = " + "{0:.2f}".format(float(badbc_count/count)) + "\n")
        f.write("goodbc_badalignment count: " + str(goodbc_badalign_count) + " = " + "{0:.2f}".format(float(goodbc_badalign_count/count)) + "\n")
        f.write("goodbc_goodalignment count: " + str(goodbc_goodalign_count) + " = " + "{0:.2f}".format(float(goodbc_goodalign_count/count)) + "\n")
        f.write("goodbc_perfectalignment count: " + str(goodbc_perfectalign_count) + " = " + "{0:.2f}".format(float(goodbc_perfectalign_count/count)) + "\n")
        
    print("Parsing of " + file + " finished at:" + str(datetime.now()) + "\n")

def parse_rna(file):
    fadapter = "GTGGTATCAACGCAGAGTACAT"
    radapter = "CTGCAGCGT"


    #alignment scoring
    match = 0
    mismatch = -1
    gap_open = -1
    gap_extend = -1
    score_dist = []

    sref_bc = set(ref_bc)
    #sref_bc_last10 = set(ref_bc_last10)
    
    exp_err_dist = []
    read_trim_len_dist = []
    prom_len_dist = []
    bc_len_dist = []
    
    count = 0
    both_adaptercount = 0
    
    lowq_count = 0
    missingadapter_count = 0
    frag_count = 0
    badbc_count = 0
    goodbc_perfectalign_count = 0
    goodbc_goodalign_count = 0
    goodbc_badalign_count = 0

    adapmismatch = defaultdict(int)
    
    lowq = defaultdict(int)
    missingadapter = defaultdict(int)
    frag = defaultdict(int)
    badbc = defaultdict(int)
    goodbc_perfectalign = defaultdict(list)
    goodbc_goodalign = defaultdict(list)
    goodbc_badalign = defaultdict(list)
    
    goodbc_perfectalign_bccounts = defaultdict(int)
    goodbc_goodalign_bccounts = defaultdict(int)
    
    print("Parsing of " + file + " started at:" + str(datetime.now()) + "\n")
    handle = gzip.open(file,"rt")
    for rec in SeqIO.parse(handle, 'fastq'):
        count +=1
        qscore = rec.letter_annotations["phred_quality"]
        read = str(rec.seq.reverse_complement())
        error = exp_err(qscore)
        exp_err_dist.append(error)
        if error <3:
            
            if fadapter in read and radapter in read:
                both_adaptercount +=1
                               
                fpos = read.find(fadapter)
                rpos = read.find(radapter)
                read_trim = read[fpos+len(fadapter):rpos]
                read_trim_len_dist.append(len(read_trim))
                #print(read_trim)
# modified
                if 15 < len(read_trim) <= 115: #figure out this line
                    prom = read_trim[4:-15] #trim off first six to get rid of the ATGGG
                    prom_len_dist.append(len(prom))
                    bc = read_trim[-12:]
                    bc_len_dist.append(len(bc))
                            
                    sbc = str(bc)
                    if sbc in sref_bc and len(prom)>0:
                        trim_prom = prom_lib[sbc][100-len(prom):]                                    
                        score = pairwise2.align.globalms(prom,trim_prom,match,mismatch,gap_open,gap_extend,score_only = True)
                        if score == 0:
                            goodbc_perfectalign_count += 1
                            tabulate(sbc,goodbc_perfectalign_bccounts)
                            tabulate_read(sbc,prom,goodbc_perfectalign)
                            
                        elif score >= (math.ceil(len(prom)*0.04)*-1):
                            goodbc_goodalign_count +=1
                            tabulate(sbc,goodbc_goodalign_bccounts)
                            tabulate_read(sbc,prom,goodbc_goodalign)
                            score_dist.append(score)
                        else:
                            score_dist.append(score)
                            goodbc_badalign_count += 1
                            tabulate_read(sbc,prom,goodbc_badalign)           
                    else:
                        badbc_count +=1
                        tabulate(read_trim,badbc)
                else:
                    frag_count +=1  
                    tabulate(read_trim,frag)
            else:
                missingadapter_count += 1
                tabulate(read,missingadapter)
            
        else:
            tabulate(read,lowq)
            lowq_count += 1
            
        #if count >50000:
            #break  
            
    combined = {}
    combined.update(goodbc_perfectalign_bccounts)
    combined.update(goodbc_goodalign_bccounts)
        
    combined_bccounts = {}
    for key in combined:
        if key in goodbc_goodalign_bccounts and key in goodbc_perfectalign_bccounts:
            combined_bccounts[key] = goodbc_goodalign_bccounts[key] + goodbc_perfectalign_bccounts[key] 
        elif key in goodbc_goodalign_bccounts:
            combined_bccounts[key] = goodbc_goodalign_bccounts[key]
        else:
            combined_bccounts[key] = goodbc_perfectalign_bccounts[key]

    print("total count: " + str(count))
    print("")
    print("lowq count: " +str(lowq_count) + " = " + "{0:.2f}".format(float(lowq_count/count)))
    print("missing_adapter count: "+ str(missingadapter_count) + " = " + "{0:.2f}".format(float(missingadapter_count/count)))
    print("frag count: " + str(frag_count) + " = " + "{0:.2f}".format(float(frag_count/count)))
    print("badbc count: " + str(badbc_count) + " = " + "{0:.2f}".format(float(badbc_count/count)))
    print("goodbc_badalignment count: " + str(goodbc_badalign_count) + " = " + "{0:.2f}".format(float(goodbc_badalign_count/count)))
    print("goodbc_goodalignment count: " + str(goodbc_goodalign_count) + " = " + "{0:.2f}".format(float(goodbc_goodalign_count/count)))
    print("goodbc_perfectalignment count: " + str(goodbc_perfectalign_count) + " = " + "{0:.2f}".format(float(goodbc_perfectalign_count/count)))
    print("")

    filename = file.split("/")[0] + "/" + file.split("/")[-1].split(".")[0]
    runname = file.split("/")[0]
    fname = file.split("/")[-1].split(".")[0]
    d = pd.DataFrame(list(combined_bccounts.items()), columns = ['Barcode', 'Counts'])
    d = d.set_index('Barcode')
    d.to_csv("05_bc_counts/" + fname + ".csv")
    d2 = pd.read_csv("05_bc_counts/" + fname + ".csv")

    BC = ""
    BC_lst = []
    for BC_rev in d2['Barcode']:
        BC = Seq(BC_rev).reverse_complement()
        BC = ''.join(BC)
        BC_lst.append(BC)
    d2['Barcode'] = BC_lst
    d2.to_csv("05_bc_counts/" + fname + ".csv", index=False)
    
# modified
    #write to files
    write_dict(combined_bccounts,"0_bccounts/" + fname + ".txt")
    write_dict(lowq,"1_lowq/" + fname + ".txt")
    write_dict(missingadapter,"2_missingadapter/" + fname + ".txt")
    write_dict(frag,"3_frag/" + fname + ".txt")
    write_dict(badbc,"4_badbc/" + fname + ".txt")
    write_dict(goodbc_badalign,"5_goodbc_badalign/" + fname + ".txt")
    write_dict(goodbc_goodalign,"6_goodbc_goodalign/" + fname + ".txt")
    write_dict(goodbc_perfectalign,"7_goodbc_perfectalign/" + fname + ".txt")

    write_dict(goodbc_goodalign_bccounts,"8_goodbc_goodalign_bccounts/" + fname + ".txt")
    write_dict(goodbc_perfectalign_bccounts,"9_goodbc_perfectalign_bccounts/" + fname + ".txt")
    
    with open("10_log_files/" + fname + ".txt","w+") as f:
        f.write("total count: " + str(count) + "\n")
        f.write("" + "\n")
        f.write("lowq count: " +str(lowq_count) + " = " + "{0:.2f}".format(float(lowq_count/count)) + "\n")
        f.write("missing_adapter count: "+ str(missingadapter_count) + " = " + "{0:.2f}".format(float(missingadapter_count/count)) + "\n")
        f.write("frag count: " + str(frag_count) + " = " + "{0:.2f}".format(float(frag_count/count)) + "\n")
        f.write("badbc count: " + str(badbc_count) + " = " + "{0:.2f}".format(float(badbc_count/count)) + "\n")
        f.write("goodbc_badalignment count: " + str(goodbc_badalign_count) + " = " + "{0:.2f}".format(float(goodbc_badalign_count/count)) + "\n")
        f.write("goodbc_goodalignment count: " + str(goodbc_goodalign_count) + " = " + "{0:.2f}".format(float(goodbc_goodalign_count/count)) + "\n")
        f.write("goodbc_perfectalignment count: " + str(goodbc_perfectalign_count) + " = " + "{0:.2f}".format(float(goodbc_perfectalign_count/count)) + "\n")

       
    print("Parsing of " + file + " finished at:" + str(datetime.now()) + "\n")


for file in DNAfilepaths:
    parse_dna(file)
    
for file in RNAfilepaths:
    parse_rna(file)