import os

di = "20200320_pJP50regall/"
rawdi = di+"raw/"
mergdi = di+"merged/"
unmergdi = di+"unmerged/"


with open("merge.sh","w+") as wfile:
	wfile.write("#!/bin/bash \n")
	for file in [f for f in os.listdir(rawdi) if os.path.isfile(rawdi + f)]:
		if "R1" in file:
			R2_file = file.replace("R1","R2")
			mergename = file.split("_")[0] + file.split("_")[1] + file.split("_")[2] + file.split("_")[3] + ".fastq.gz"
			#mergename = file.split("_")[0] + ".fastq.gz"
			#mergename = file.split("_")[0] + file.split("_")[1] + file.split("_")[2] + ".fastq.gz"
			#if "DNA_I" in file:
			#	mergename = file.split("_")[0] + file.split("_")[1] + file.split("_")[2] + ".fastq.gz"
			#if "DNA_V" in file:
			#	mergename = file.split("_")[0] + file.split("_")[1] + file.split("_")[2] + file.split("_")[3] + ".fastq.gz"
			#if "RNA_I" in file:
			#	mergename = file.split("_")[0] + file.split("_")[1] + file.split("_")[2] + file.split("_")[3] + file.split("_")[4] + ".fastq.gz"
			#if "RNA_V" in file:
			#	mergename = file.split("_")[0] + file.split("_")[1] + file.split("_")[2] + file.split("_")[3] + file.split("_")[4] + file.split("_")[5] + file.split("_")[6] + ".fastq.gz"
			#if "VIVO" in file:
			#	mergename = file.split("_")[0] + file.split("_")[1] + file.split("_")[2] + file.split("_")[3] + ".fastq.gz"


			R1_unmerged = file.replace(".fastq.gz","_unmerged.fastq.gz")
			R2_unmerged = R2_file.replace(".fastq.gz","_unmerged.fastq.gz")
			wfile.write("source ./bbmap/bbmerge.sh in1=" + rawdi + file + " in2=" +rawdi + R2_file +
			" out=" + mergdi + mergename + " outu1=" + unmergdi + R1_unmerged + " outu2=" + unmergdi + R2_unmerged + "\n")


