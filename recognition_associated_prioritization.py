import getopt,sys,subprocess,os
import csv
import multiprocessing
import pandas as pd
import numpy as np
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.SeqIO.FastaIO import SimpleFastaParser
from math import log, exp

def getR(neo_seq,iedb_seq): 
    align_score = []
    a = 26
    k = 4.86936
    for seq in iedb_seq:
        aln_score = aligner(neo_seq,seq)
        # print(aln_score)
        if aln_score:
            localds_core = max([line[2] for line in aln_score])
            align_score.append(localds_core)

    bindingEnergies = list(map(lambda x: -k * (a - x), align_score))
    lZk = logSum(bindingEnergies + [0])
    lGb = logSum(bindingEnergies)
    R=exp(lGb-lZk)
    return R

def getiedbseq(iedb_path):
    iedb_seq = []
    with open(iedb_path, 'r') as fin: #'/data8t_2/zzt/antigen.garnish/iedb.fasta'
        for t, seq in SimpleFastaParser(fin):
            iedb_seq.append(seq)
    return iedb_seq

def logSum(v):
    ma = max(v)
    return log(sum(map(lambda x: exp(x-ma),v))) + ma

def aligner(seq1,seq2):
    matrix = matlist.blosum62
    gap_open = -11
    gap_extend = -1
    aln = pairwise2.align.localds(seq1.upper(), seq2.strip().split('+')[0].upper(), matrix, gap_open, gap_extend)
    return aln
    
def write_file(a_list, name):
    textfile = open(name, "w")
    for element in a_list:
        textfile.write(element + "\n")
    textfile.close()
    return

def get_wt_bindaff(wt_seq,hla,output_folder,netmhc_path):
    # write_file(wt_seq, "stage.pep")
    os.system("mkdir "+output_folder+"/tmp")
    with open(output_folder+"/tmp/wt.pep", "w") as pepfile:
        pepfile.write("%s"% wt_seq)
    pepfile.close()
    args = netmhc_path+" -p "+output_folder+"/tmp/wt.pep -a "+ hla+" -l "+str(len(wt_seq))+" -BA >> "+output_folder+"/tmp/wt.csv"
    subprocess.call(args, shell=True)  
    wt_bindaff = 1
    with open(output_folder+"/tmp/wt.csv") as f:
        data = f.read()
    nw_data = data.split('-----------------------------------------------------------------------------------\n')
    for i in range(len(nw_data)):
        if i%4 == 2:
            wt_bindaff = nw_data[i].strip().split()[15]
            break
    os.system("rm -rf "+output_folder+"/tmp")
    return wt_bindaff


def main():
    opts,args=getopt.getopt(sys.argv[1:],"hi:I:o:n:f:a:t:p:",["input_folder=","iedb_fasta=","output_folder=","netmhc_path=","foreignness_score=","agretopicity=","alteration_type=", "prefix="])
    input_folder=""
    iedb_fasta=""
    output_folder=""
    netmhc_path=""
    foreignness_score= 1e-16
    agretopicity=0.1
    alteration_type="snv,indel,fusion,splicing"
    prefix=""
    USAGE='''
        This script rank neoantigens by recognition associated evaluation (foreigness and agretopicity)
        usage: python recognition_associated_prioritization.py -i <input_folder> -o <output_folder> \
            -f <foreignness_score> -a <agretopicity> -t <alteration_type> -p <prefix>
            required argument:
                -i | --input_folder : input folder including result file from bindstab output
                -I | --iedb_fasta : path to iedb fasta reference file
                -o | --output_folder : output folder to store result
                -n | --netmhc_path : path to run netmhcpan
                -f | --foreignness_score : neoantigen threshold for neoantigen
                -a | --agretopicity : neoantigen threshold for agretopicity
                -t | --alteration_type : neoantigen from alteration type to rank (default is "snv,indel,fusion,splicing")
                -p | --prefix : prefix of output file
    '''
    for opt,value in opts:
        if opt =="h":
            print (USAGE)
            sys.exit(2)
        elif opt in ("-i","--input_folder"):
            input_folder=value
        elif opt in ("-I","--iedb_fasta"):
            iedb_fasta=value
        elif opt in ("-o","--output_folder"):
            output_folder =value 
        elif opt in ("-n","--netmhc_path"):
            netmhc_path =value 
        elif opt in ("-f","--foreignness_score"):
            foreignness_score =float(value) 
        elif opt in ("-a","--agretopicity"):
            agretopicity =float(value) 
        elif opt in ("-t","--alteration_type"):
            alteration_type =value 
        elif opt in ("-p","--prefix"):
            prefix =value 

    if (input_folder =="" or iedb_fasta=="" or output_folder=="" or netmhc_path==""):
        print (USAGE)
        sys.exit(2)		

    wt_bindaff_list = []
    wt_list = []
    tmp_wt_bindaff_file =csv.reader(open(input_folder+"/tmp_identity/"+prefix+"_bindaff_filtered.tsv"), delimiter="\t")
    for line in tmp_wt_bindaff_file:
        if line[7] != "":
            wt_bindaff_list.append(line[7])
            wt_list.append(line[2])
    wt_bindaff_list.pop(0)
    wt_list.pop(0)

    snv_indel_file = open(output_folder+"/../info/"+prefix+"_snv_indel.vcf")
    if os.path.exists(output_folder+"/../info/"+prefix+"_fusion.tsv"):
        fusion_file = open(output_folder+"/../info/"+prefix+"_fusion.tsv")
    else:
        fusion_file = []
    if os.path.exists(output_folder+"/../info/"+prefix+"_splicing.csv"):
        splicing_file = open(output_folder+"/../info/"+prefix+"_splicing.csv")
    else:
        splicing_file = []
    snv_indel = []
    fusion = []
    splicing = []

    for line in snv_indel_file:
        snv_indel.append(str(line))
    for line in fusion_file:
        fusion.append(str(line))
    for line in splicing_file:
        splicing.append(str(line))


    iedb_seq = getiedbseq(iedb_fasta)
    reader = csv.reader(open(input_folder+"/"+prefix+"_candidate_pmhc.csv"), delimiter=",")
    fields=next(reader)
    fields.append("Foreignness")
    fields.append("Agretopicity")
    fields.append("SourceAlterationDetail")
    data_raw = []
    data_exist = [] # save existing hla, mutant_type peptide
    agre_exist = []
    for line in reader:
        R = getR(line[1], iedb_seq)
        line.append(R)
        mt_bindaff = float(line[3])
        identity = line[5]
        if "SNV" in identity:
            while (wt_list.pop(0)!=line[2]):
                wt_bindaff_list.pop(0)
            wt_bindaff = float(wt_bindaff_list.pop(0))
        else:
            wt_bindaff = 0.001 # set to a small number for large agretopicity (only rank by foreigness)
        
        A = mt_bindaff/float(wt_bindaff) 
        if ([line[0],line[1]] in data_exist):
            index = data_exist.index([line[0],line[1]])
            if (A < agre_exist[index]):
                agre_exist[index] = A
                data_raw[index][8]=A
        else:
            data_exist.append([line[0],line[1]])
            agre_exist.append(A)
            line.append(A)
            data_raw.append(line)

        line_info_string = ""
        if (identity.strip().split('_')[0]=="SNV" or identity.strip().split('_')[0]=="INDEL"):
            line_num = int(identity.strip().split('_')[1])
            snv_indel_line = snv_indel[line_num-1]
            ele = snv_indel_line.strip().split('\t')
            if len(ele) == 14: # annotation software is vep
                annotation_info = ["Uploaded_variation","Location","Allele","Gene","Feature","Feature_type",
                                    "Consequence","cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","Extra"]
                for i in range(0,len(ele),1):
                    line_info_string+=annotation_info[i]+"$"+ele[i]+"#"
            elif len(ele)==11:
                annotation_info = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","normal","tumor"]
                for i in range(0,len(ele),1):
                    line_info_string+=annotation_info[i]+"$"+ele[i]+"#"
            else:
                continue
        elif (identity.strip().split('_')[0]=="FUSION"):
            line_num = int(identity.strip().split('_')[1])
            fusion_line = fusion[line_num-1]
            ele = fusion_line.strip().split('\t')
            annotation_info = ["FusionName","JunctionReadCount","SpanningFragCount","est_J","est_S","SpliceType","LeftGene","LeftBreakpoint",
                                "RightGene","RightBreakpoint","LargeAnchorSupport","FFPM","LeftBreakDinuc","LeftBreakEntropy","RightBreakDinuc",
                                "RightBreakEntropy","annots","CDS_LEFT_ID","CDS_LEFT_RANGE","CDS_RIGHT_ID","CDS_RIGHT_RANGE","PROT_FUSION_TYPE",
                                "FUSION_MODEL","FUSION_CDS","FUSION_TRANSL","PFAM_LEFT","PFAM_RIGHT"]
            for i in range(0, len(ele),1):
                line_info_string+=annotation_info[i]+"$"+ele[i]+"#"
        elif (identity.strip().split('_')[0]=="SP"):
            line_num = int(identity.strip().split('_')[1])
            splicing_line = splicing[line_num-1]
            ele = splicing_line.strip().split('\t')
            annotation_info = ["chrom","txStart","txEnd","isoform","protein","strand","cdsStart","cdsEnd","gene","exonNum",
                                "exonLens","exonStarts","ensembl_transcript"]
            for i in range(0, len(ele),1):
                line_info_string+=annotation_info[i]+"$"+ele[i]+"#"
        else:
            continue
        line[5] = identity.strip().split('_')[0]
        line.append(line_info_string)

    picked_rows = []
    alt_type = alteration_type.replace(" ", "").strip().split(",")
    for line in data_raw:
        type = line[5].strip().split("_")[0]
        if (type=="SP"):
            type="SPLICING"
        if type.lower() in alt_type:
            picked_rows.append(line)
    data=pd.DataFrame(picked_rows)
    data.columns=fields
    data.BindAff = data.BindAff.astype(float)
    data.Foreignness = data.Foreignness.astype(float)
    data.Agretopicity = data.Agretopicity.astype(float)
    data_filtered_true = data[(data.Foreignness>foreignness_score) | (data.Agretopicity<agretopicity)]#(data.filtered==True)
    data_filtered_false = data[~((data.Foreignness>foreignness_score) | (data.Agretopicity<agretopicity))]
    data_filtered_true["rank_val"]=data_filtered_true["BindAff"].rank(pct = True)
    data_filtered_false["rank_val"]=data_filtered_false["BindAff"].rank(pct = True)+1
    data=pd.concat([data_filtered_true,data_filtered_false])
    data["Rank"]=data["rank_val"].rank(method="first")
    data=data.sort_values("Rank")
    data=data.astype({"Rank":int})
    data=data.drop(['BindLevel', 'rank_val'], axis=1)
    data.loc[data['Identity'] != "SNV","WT_pep"] = "-"
    data.loc[data['Identity'] != "SNV","Agretopicity"] = "-"
    data.to_csv(output_folder+"/"+prefix+"_neoantigen_rank_recognition_associated.tsv",header=1,sep='\t',index=0)

if __name__ == '__main__':
    main()