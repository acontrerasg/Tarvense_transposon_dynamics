#!/usr/bin/python3
## Author: Dario Galanti, Dec 2020
## Input: simil unionbed file with methylation info of TIRs flanks.
## Input Headers: chr st end flk_chr flk_st flk_end TIP_ID fam flk(1/2) n_spls_with spls_with spl1 spl2 spl3…
## Aim: Calculate coordinates relative to flank and methylation of samples with and without TIP
## Output: file with methylation info of TIP flanks of samples with vs without insertion.
## Output headers: pos flk(1/2) met_splsWith met_splsWithout
## Run: python3 TIPflks_met.py input.bed output.txt
## Run local: python TIPflks_met.py TA_R_unionTIPflanks.bed TA_R_TIPflanks_met.txt

## Dependencies: conda activate Python

## Process:
## 1) Retrieve IDs of samples with and without the insertion for each TIP

## Import modules
import sys

## Define input and create output directory
fin = str(sys.argv[1])
fout = str(sys.argv[2])
MSC = 2		## Minor State Count: Minimum number of samples with different state (insertion/no insertion)


## 1) Retrieve IDs of samples with and without the insertion for each TIP
line_num = 0
with open(fin, "r") as fin, open(fout, "w") as fout:
	for line in fin:
		line_num += 1
		line = line.strip()
		line = line.split("\t")
		if line_num == 1:
			spls = line[11:]
			print("pos","ID","fam","flk","met_with","met_without",sep="\t", file=fout)		# Print headers
		else:
		# Retrieve spls with insertion/deletion
			spls_with = line[10].split(",")
		# Calulate methylation of samples with and without
			spl_index = 10
			n_spls_with = 0
			n_spls_without = 0
			met_with = 0
			met_without = 0
			for s in spls:
				spl_index += 1
				if line[spl_index] != "NA":
					if s in spls_with:
						n_spls_with += 1
						met_with += float(line[spl_index])
					else:
						n_spls_without += 1
						met_without += float(line[spl_index])
		# If there is data for both samples with and without insertion/deletion, continue
			if n_spls_with > MSC and n_spls_without > MSC:
				met_with = round((met_with/n_spls_with),2)
				met_without = round((met_without/n_spls_without),2)
		# Calculate new coordinates (relative to the flank)
				if line[8] == "flk1":
					pos = int(line[2]) - int(line[5])	# Position - bin end
				else:
					pos = int(line[2]) - int(line[4])	# Position - bin start
		# Print line (if not NA)
				print(pos,line[6],line[7],line[8],met_with,met_without,sep="\t", file=fout)
				


