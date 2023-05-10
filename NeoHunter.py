import os,sys
import multiprocessing
import subprocess
import yaml
import warnings
import argparse
import csv
import time
import logging
import pandas as pd
import psutil
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
warnings.filterwarnings("ignore")

def output_time(string, output_folder):
    textfile = open(output_folder+"/../TimeStamp.txt", "a")
    textfile.write(string + "\n")
    textfile.close()
    return

def filter_vcf(prefix, input_vcf_file, output_folder, tumor_depth_threshold, tumor_vaf_threshold, normal_vaf_threshold):
    f_filter=open(output_folder+"/"+prefix+"_filter.vcf",'w')
    input_file_path = os.path.join(output_folder,input_vcf_file)
    for line in open(input_file_path):
        if line.startswith("#"):
            f_filter.write(line)
        else:
            record=line.strip().split('\t')
            tumor_info=record[10]
            normal_info=record[9]
            tumor_depth=tumor_info.split(':')[1].split(',')[1]
            tumor_vaf=tumor_info.split(':')[2].split(',')
            normal_vaf=normal_info.split(':')[2].split(',')
            for i in range(0,len(tumor_vaf),1):
                if int(tumor_depth)>=int(tumor_depth_threshold) and float(tumor_vaf[i])>=float(tumor_vaf_threshold) and float(normal_vaf[i])<=float(normal_vaf_threshold):
                    f_filter.write(line)

    f_filter.close()


def fastq2bam(prefix, fastq_1_path, fastq_2_path, output_folder, reference, bwa_path, thread):
    cmd_bwa=bwa_path + " mem -t "+str(thread)+" " + reference + " " + fastq_1_path + " " +fastq_2_path + " > " + output_folder+"/tmp_"+prefix+".sam"
    logging.info(cmd_bwa)
    subprocess.call(cmd_bwa, shell=True, executable="/bin/bash")


def alignment(prefix, bam_path, output_folder, log_folder, samtools_path, thread):
    cmd_samtools_view=samtools_path + " view -bhS -@ "+str(thread)+" " + bam_path +" -o " + output_folder+"/"+"tmp_"+ prefix+".bam > " + log_folder + "/" + prefix + "_samtools_view.log 2>&1"
    logging.info(cmd_samtools_view)
    subprocess.call(cmd_samtools_view, shell=True, executable="/bin/bash")
    cmd_samtools_sort=samtools_path + " sort -@ "+str(thread)+" " + output_folder+"/"+"tmp_"+ prefix +".bam -o " + output_folder+"/"+ prefix+".bam > " + log_folder + "/" + prefix + "_samtools_sort.log 2>&1"
    logging.info(cmd_samtools_sort)
    subprocess.call(cmd_samtools_sort, shell=True, executable="/bin/bash")
    cmd_samtools_index=samtools_path + " index " + output_folder+"/"+ prefix+".bam > " + log_folder + "/" + prefix + "_samtools_index.log 2>&1"
    logging.info(cmd_samtools_index)
    subprocess.call(cmd_samtools_index, shell=True, executable="/bin/bash")
    cmd_samtools_flagstat =samtools_path + " flagstat -@ "+str(thread)+" " + output_folder+"/"+ prefix+".bam > "+output_folder+"/"+prefix+"_stat.tsv"
    logging.info(cmd_samtools_flagstat)
    subprocess.call(cmd_samtools_flagstat, shell=True, executable="/bin/bash")

def pre_mutect2(prefix, reference, input_folder, output_folder, log_folder, OneKG_path, mills_path, dbsnp_path, gatk_path, picard_path, samtools_path):
    cmd_picard="java -Xmx4G -jar " + picard_path + " MarkDuplicates INPUT=" + input_folder+"/"+ prefix+".bam OUTPUT=" + output_folder+"/"+ prefix+"_mkdup_filter.bam METRICS_FILE=" + output_folder+"/"+prefix+"_dup_qc.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT > " + log_folder + "/" + prefix + "_markdup.log 2>&1"
    logging.info(cmd_picard)
    subprocess.call(cmd_picard, shell=True, executable="/bin/bash")
    cmd_samtools_index=samtools_path + " index " + output_folder+"/"+ prefix+"_mkdup_filter.bam"
    logging.info(cmd_samtools_index)
    subprocess.call(cmd_samtools_index, shell=True, executable="/bin/bash")
    cmd_add_readgroup="java -Xmx4G -jar " + picard_path + " AddOrReplaceReadGroups I=" + output_folder+"/"+ prefix+"_mkdup_filter.bam O=" + output_folder+"/"+ prefix+"_mkdup_filter_add.bam SO=coordinate VALIDATION_STRINGENCY=SILENT RGID=" +  prefix +  " RGLB=" +  prefix + " RGPL=illumina RGSM="+  prefix + " RGPU=NextSeq > " + log_folder + "/" + prefix + "_addreadgroup.log 2>&1"
    logging.info(cmd_add_readgroup)
    subprocess.call(cmd_add_readgroup, shell=True, executable="/bin/bash")
    cmd_buildbamindex="java -Xmx4G -jar " + picard_path + " BuildBamIndex I=" + output_folder+"/"+ prefix+"_mkdup_filter_add.bam O=" + output_folder+"/"+ prefix+"_mkdup_filter_add.bam.bai VALIDATION_STRINGENCY=SILENT > " + log_folder + "/" +  prefix + "_buildindex.log 2>&1"
    logging.info(cmd_buildbamindex)
    subprocess.call(cmd_buildbamindex, shell=True, executable="/bin/bash")
    cmd_BaseRecalibrator="java -Xmx4G -jar " + gatk_path + " BaseRecalibrator -R " + reference + " -I " + output_folder+"/"+ prefix+"_mkdup_filter_add.bam --known-sites " + OneKG_path + " --known-sites " + dbsnp_path + " -O " + output_folder + "/" + prefix + ".table > " + log_folder + "/" +  prefix + "_BaseRecalibrator.log 2>&1"
    logging.info(cmd_BaseRecalibrator)
    subprocess.call(cmd_BaseRecalibrator, shell=True, executable="/bin/bash")
    cmd_PrintReads="java -Xmx4G -jar " + gatk_path + " ApplyBQSR -R " + reference + " -I " + output_folder+"/"+ prefix+"_mkdup_filter_add.bam --bqsr-recal-file " + output_folder + "/" + prefix + ".table -O " + output_folder + "/" + prefix + "_recal.bam > " + log_folder + "/" +  prefix + "_ApplyBQSR.log 2>&1"
    logging.info(cmd_PrintReads)
    subprocess.call(cmd_PrintReads, shell=True, executable="/bin/bash")


def mutation_detection(prefix, reference, output_folder, log_folder, gatk_path, interval_list_path,  tumor_depth, tumor_vaf, normal_vaf):
    if interval_list_path != "":
        cmd_gatk= "java -jar "+ gatk_path + " Mutect2 \
            -R "+reference+" \
            -I "+output_folder+"/"+prefix+"_tumor_recal.bam \
            -I "+output_folder+"/"+prefix+"_normal_recal.bam \
            -O "+output_folder+"/"+prefix+"_mutect2.vcf \
            -L "+ interval_list_path+" \
            > "+log_folder+"/"+prefix+"_mutect2.log 2>&1"
        logging.info(cmd_gatk)
        subprocess.call(cmd_gatk, shell=True, executable="/bin/bash")
    else :
        cmd_gatk= "java -jar "+ gatk_path + " Mutect2 \
            -R "+reference+" \
            -I "+output_folder+"/"+prefix+"_tumor_recal.bam \
            -I "+output_folder+"/"+prefix+"_normal_recal.bam \
            -O "+output_folder+"/"+prefix+"_mutect2.vcf \
            > "+log_folder+"/"+prefix+"_mutect2.log 2>&1"
        logging.info(cmd_gatk)
        subprocess.call(cmd_gatk, shell=True, executable="/bin/bash")
    filter_vcf(prefix, prefix+"_mutect2.vcf", output_folder, tumor_depth, tumor_vaf, normal_vaf)


def hla_typing(prefix, rna_fastq_path_1, rna_fastq_path_2, output_folder, log_folder, optitype_path):
    start_t = time.time()
    logging.info("Start HLA Typing")
    hla_output_path = os.path.join(output_folder, "hla_typing")
    # run optitype for hla typings
    if not os.path.exists(hla_output_path):
        os.mkdir(hla_output_path)
    cmd_hla = "conda run -n Optitype_env python " + optitype_path + " -i " + rna_fastq_path_1 + " "+rna_fastq_path_2+" --rna -o " + hla_output_path+" > "+log_folder+"/"+prefix+"_hla_typing.log 2>&1"
    logging.info(cmd_hla)
    subprocess.call(cmd_hla, shell=True, executable="/bin/bash")
    logging.info("Finish HLA Typing")
    # move the directory
    tmp_dir=os.listdir(hla_output_path)
    hla_tsv = hla_output_path+"/"+tmp_dir[0]+"/"+tmp_dir[0]+"_result.tsv"
    hla_pdf = hla_output_path+"/"+tmp_dir[0]+"/"+tmp_dir[0]+"_coverage_plot.pdf"
    cmd_rename="mv "+hla_tsv+" "+hla_output_path+"/"+prefix+"_hla.tsv"
    subprocess.call(cmd_rename, shell=True, executable="/bin/bash")
    cmd_rename="mv "+hla_pdf+" "+hla_output_path+"/"+prefix+"_hla_coverage_plot.pdf"
    subprocess.call(cmd_rename, shell=True, executable="/bin/bash")
    cmdrm = "rm -rf "+hla_output_path+"/"+tmp_dir[0]
    subprocess.call(cmdrm, shell=True, executable="/bin/bash")
    end_t = time.time()
    output_time("Time Takes for HLA typing is "+str(round(end_t-start_t,2))+" s", hla_output_path)


def kallisto_expression(rna_fastq_1,rna_fastq_2,kallisto_cdna_path,kallisto_path,output_folder,log_folder,thread):
    start_t = time.time()
    cdna_path_dir = os.path.dirname(kallisto_cdna_path)
    cnd_file_prefix = os.path.splitext(os.path.basename(kallisto_cdna_path))[0]
    kallisto_index_path = cdna_path_dir + "/" + cnd_file_prefix + ".idx"
    cmd_kallisto_index = kallisto_path + " index -i " + kallisto_index_path + " " + kallisto_cdna_path + " > " +  log_folder + "/kallisto_index.log 2>&1"
    logging.info(cmd_kallisto_index)
    subprocess.call(cmd_kallisto_index, shell=True, executable="/bin/bash")
    if os.path.exists(rna_fastq_2):
        cmd_kallisto_quant = kallisto_path + " quant -i " + kallisto_index_path + " -t "+str(thread)+" -b 100 -o " + output_folder + " " + rna_fastq_1 + " " + rna_fastq_2 + " > " +  log_folder + "/kallisto.log 2>&1"
        logging.info(cmd_kallisto_quant)
        subprocess.call(cmd_kallisto_quant, shell=True, executable="/bin/bash")
    else:
        cmd_kallisto_quant = kallisto_path + " quant -i " + kallisto_index_path + " -t "+str(thread)+" -b 100 --single -l 200 -s 30 -o " + output_folder + " " + rna_fastq_1 +  " > " +  log_folder + "/kallisto.log 2>&1"
        logging.info(cmd_kallisto_quant)
        subprocess.call(cmd_kallisto_quant, shell=True, executable="/bin/bash")
    end_t = time.time()
    output_time("Time Takes for Transcript Quantification is "+str(round(end_t-start_t,2))+" s",log_folder)
    return


def snv_indel_pred(prefix, input_vcf, reference, human_peptide_path, tumor_abundance, output_folder, perl_path, gatk_path, vep_path, vep_cache, snpeff_path, funcotator_source, annotation_software):
    start_t = time.time()
    logging.info("Running SNV_Indel_prediction with annotation software: "+annotation_software)
    snv_indel_folder = os.path.join(output_folder,prefix+"_snv_indel_output")
    if (annotation_software == "VEP"):
        os.environ["PERL5LIB"] = "/home/zzt/anaconda3/envs/zzt/lib/site_perl/5.26.2/x86_64-linux-thread-multi:$PERL5LIB"
        cmd_vep= "conda run -n zzt "+perl_path+" "+vep_path + " --species homo_sapiens --assembly \
            GRCh37 --no_stats --buffer_size 5000 --sift b --ccds --uniprot --hgvs \
            --symbol --numbers --domains --gene_phenotype --canonical --protein \
            --biotype --tsl --variant_class --shift_hgvs 1 \
            --check_existing --total_length --allele_number --no_escape \
            --xref_refseq --failed 1 --flag_pick_allele --pick_order \
            canonical,tsl,biotype,rank,ccds,length --dir \
            "+vep_cache+" --fasta "+reference+" \
            --format vcf --input_file "+ input_vcf +" \
            --output_file \
            "+snv_indel_folder+"/"+prefix+"_snv_indel_annotation.vcf \
            --offline --pubmed --fork 4 --polyphen b --af --af_1kg --af_esp \
            --af_gnomad --regulatory"
        logging.info(cmd_vep)
        subprocess.call(cmd_vep, shell=True, executable="/bin/bash")
        os.environ["PERL5LIB"] = "/home/zzt/perl5/lib/perl5"
    elif (annotation_software == "SnpEff"):
        cmd_snpeff = "java -Xmx4g -jar " + snpeff_path + " GRCh37.75 \
            "+input_vcf+" > "+snv_indel_folder+"/"+prefix+"_snv_indel_annotation.vcf"
        logging.info(cmd_snpeff)
        subprocess.call(cmd_snpeff, shell=True, executable="/bin/bash")
    elif (annotation_software == "Funcotator"):
        cmd_funcotator = gatk_path+" Funcotator --variant "+input_vcf+" \
                --reference "+reference+"\
                --ref-version hg19 \
                --data-sources-path "+ funcotator_source+" \
                --output "+snv_indel_folder+"/"+prefix+"_snv_indel_annotation.vcf \
                --output-file-format VCF"
        logging.info(cmd_funcotator)
        subprocess.call(cmd_funcotator, shell=True, executable="/bin/bash")
    else:
        logging.error("Type of annotation software is invalid")
    end_t = time.time()
    transcript_quantification_folder = os.path.join(output_folder,prefix+"_transcript_quantification")
    output_time("Time Takes for SNV/INDEL Annotation is "+str(round(end_t-start_t,2))+" s", output_folder)
    cmd_mutation_peptide="python ./annotation2fasta.py -i " + snv_indel_folder + "/"+prefix+"_snv_indel_annotation.vcf \
                    -o " + output_folder + " -p " + human_peptide_path + " \
                    -r " + reference+" \
                    -s " + annotation_software + " -e " +transcript_quantification_folder+"/abundance.tsv \
                    -t "+ str(tumor_abundance) + " -P "+prefix
    logging.info(cmd_mutation_peptide)
    subprocess.call(cmd_mutation_peptide, shell=True, executable="/bin/bash")


def fusion_pred(prefix, RNA_tumor_1, RNA_tumor_2, star_fusion_dataset, output_folder, star_fusion_path, star_path, tumor_abundance, min_FFPM, thread):
    start_t = time.time()
    star_fusion_output_path = os.path.join(output_folder, prefix+"_fusion_output")
    if (RNA_tumor_2 == "" or RNA_tumor_2 == "None"):
        run_star_fusion = "conda run -n base "+star_fusion_path+" \
                --genome_lib_dir "+star_fusion_dataset+" \
                --examine_coding_effect --left_fq "+RNA_tumor_1+" \
                --output_dir "+star_fusion_output_path+" \
                --STAR_PATH "+star_path+" \
                --min_FFPM "+str(min_FFPM)+" \
                --CPU "+str(thread)
        logging.info(run_star_fusion)
        subprocess.call(run_star_fusion, shell=True, executable="/bin/bash")
    else:
        run_star_fusion = "conda run -n base "+star_fusion_path+" \
                --genome_lib_dir "+star_fusion_dataset+" \
                --examine_coding_effect --left_fq "+RNA_tumor_1+" \
                --right_fq "+RNA_tumor_2+" \
                --output_dir "+star_fusion_output_path+" \
                --STAR_PATH "+star_path+" \
                --min_FFPM "+str(min_FFPM)+" \
                --CPU "+str(thread)
        logging.info(run_star_fusion)
        subprocess.call(run_star_fusion, shell=True, executable="/bin/bash")
    transcript_quantification_folder = os.path.join(output_folder,prefix+"_transcript_quantification")
    parse_fusion_pred = "python parse_star_fusion.py \
        -i "+star_fusion_output_path+"/star-fusion.fusion_predictions.abridged.coding_effect.tsv \
        -e "+transcript_quantification_folder+"/abundance.tsv \
        -o "+output_folder+" \
        -p "+prefix+" \
        -t "+str(tumor_abundance)
    logging.info(parse_fusion_pred)
    subprocess.call(parse_fusion_pred, shell=True, executable="/bin/bash")
    end_t = time.time()
    output_time("Time Takes for Gene Fusion Detection is "+str(round(end_t-start_t,2))+" s", output_folder)


def splicing_pred(prefix, RNA_tumor_1, RNA_tumor_2, refseq_ann, star_path, star_fusion_dataset, asneo_path,output_folder, asneo_ref_genome, tumor_abundance, thread):
    start_t = time.time()
    splicing_output_path = os.path.join(output_folder, prefix+"_splicing_output/")
    rna_alignment_path = os.path.join(splicing_output_path, prefix+"_rna_alignment_output/")
    if (RNA_tumor_1.endswith('.gz')):
        cmd_star = "conda run -n base "+star_path+" --genomeDir "+star_fusion_dataset+"/ref_genome.fa.star.idx \
                    --readFilesIn "+RNA_tumor_1+" "+RNA_tumor_2+ " --runThreadN "+str(thread)+" –outFilterMultimapScoreRange 1 \
                    --outFilterMultimapNmax 20 –outFilterMismatchNmax 10 --alignIntronMax 500000 \
                    –alignMatesGapMax 1000000 --sjdbScore 2 –alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory \
                    –outFilterMatchNminOverLread 0.33 –outFilterScoreMinOverLread 0.33 --sjdbOverhang 150 \
                    --outSAMstrandField intronMotif –sjdbGTFfile "+refseq_ann +" --readFilesCommand 'gunzip -c' \
                    --outFileNamePrefix "+rna_alignment_path
        logging.info(cmd_star)
        subprocess.call(cmd_star, shell=True, executable="/bin/bash")
    else:
        cmd_star = "conda run -n base "+star_path+" --genomeDir "+star_fusion_dataset+"/ref_genome.fa.star.idx \
                    --readFilesIn "+RNA_tumor_1+" "+RNA_tumor_2+ " --runThreadN "+str(thread)+" –outFilterMultimapScoreRange 1 \
                    --outFilterMultimapNmax 20 –outFilterMismatchNmax 10 --alignIntronMax 500000 \
                    –alignMatesGapMax 1000000 --sjdbScore 2 –alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory \
                    –outFilterMatchNminOverLread 0.33 –outFilterScoreMinOverLread 0.33 --sjdbOverhang 150 \
                    --outSAMstrandField intronMotif –sjdbGTFfile "+refseq_ann +" \
                    --outFileNamePrefix "+rna_alignment_path
        logging.info(cmd_star)
        subprocess.call(cmd_star, shell=True, executable="/bin/bash")
    transcript_quantification_folder = os.path.join(output_folder,prefix+"_transcript_quantification")
    tmp_fasta_folder = os.path.join(output_folder,prefix+"_tmp_fasta")
    cmd_asneo = "conda run -n ASNEO_env python "+asneo_path+"/ASNEO.py -j "+rna_alignment_path+"/SJ.out.tab \
                -g "+asneo_ref_genome+" -o "+tmp_fasta_folder+" -l 8,9,10,11 -p "+prefix +" -t "+ str(tumor_abundance)+" \
                -e "+transcript_quantification_folder+"/abundance.tsv"
    logging.info(cmd_asneo)
    subprocess.call(cmd_asneo, shell=True, executable="/bin/bash")
    end_t = time.time()
    cmd_cat = "cat "+tmp_fasta_folder+"/"+prefix+"_splicing_*  > "+output_folder+"/"+prefix+"_splicing.fasta"
    logging.info(cmd_cat)
    os.system(cmd_cat)
    output_time("Time Takes for Splicing Variant Detection is "+str(round(end_t-start_t,2))+" s", output_folder)


def snv_indel_fusion_netmhc(prefix, netmhc_path, input_folder, output_folder, hla_str):
    run_netmhc = netmhc_path+" -a "+hla_str+" -f "+input_folder+"/"+prefix+"_snv_indel_fusion.fasta -l '8,9,10,11' -BA > "+output_folder+"/"+hla_str+"_all_tmp_hla_netmhc.txt"
    logging.info(run_netmhc)
    subprocess.call(run_netmhc, shell=True, executable="/bin/bash")


def splicing_netmhc(prefix, netmhc_path, input_folder, output_folder, hla_str, length):
    run_netmhc = netmhc_path+" -a "+hla_str+" -f "+input_folder+"/"+prefix+"_splicing_"+length+".fasta -l "+length+" -BA > "+output_folder+"/"+hla_str+"_all_tmp_hla_netmhc"+length+".txt"
    logging.info(run_netmhc)
    subprocess.call(run_netmhc, shell=True, executable="/bin/bash")


def snv_indel_wt_netmhc(prefix, netmhc_path, input_folder, output_folder, hla):
    run_netmhc = netmhc_path+" -a "+hla+" -f "+input_folder+"/"+prefix+"_snv_indel_wt.fasta -l '8,9,10,11' -BA > "+output_folder+"/tmp_netmhc_wt/"+hla+"_wt_tmp_hla_netmhc.txt"
    print(run_netmhc)
    subprocess.call(run_netmhc, shell=True, executable="/bin/bash")

def mutation_netmhc_parallel(prefix, netmhc_path, input_folder, output_folder, hla_str):
    hla_list = list(hla_str.strip().split(","))
    logging.info(hla_list)
    target_dir = output_folder+"/tmp_netmhc"
    if not os.path.exists(target_dir):
        os.mkdir(target_dir)
    else:
        os.removedirs(target_dir)
        os.mkdir(target_dir)
    os.mkdir(output_folder+"/tmp_netmhc_wt")
    netmhc_hla_process=[]
    for hla in hla_list:
        run_netmhc = multiprocessing.Process(target=snv_indel_fusion_netmhc,args=(prefix, netmhc_path, input_folder, target_dir, hla))#"./netmhc_parallel.sh "+netmhc_path+" "+output_folder+" "+hla
        netmhc_hla_process.append(run_netmhc)
        for length in [8,9,10,11]:
            run_splicing_netmhc = multiprocessing.Process(target=splicing_netmhc,args=(prefix, netmhc_path, input_folder, target_dir, hla, str(length)))
            netmhc_hla_process.append(run_splicing_netmhc)
        run_netmhc_wt = multiprocessing.Process(target=snv_indel_wt_netmhc,args=(prefix, netmhc_path, input_folder, output_folder, hla))#"./netmhc_parallel.sh "+netmhc_path+" "+output_folder+" "+hla
        netmhc_hla_process.append(run_netmhc_wt)
    for p in netmhc_hla_process:
        p.daemon = True
        p.start()
    for p in netmhc_hla_process:
        p.join()
    # concatenate all netmhc output file to a single file
    os.system("cat "+target_dir+"/* > "+output_folder+"/"+prefix+"_bindaff_raw.tsv")
    os.system("rm -rf "+target_dir)
    os.system("cat "+output_folder+"/tmp_netmhc_wt/* > "+output_folder+"/"+prefix+"_snv_indel_bindaff_wt.tsv")
    os.system("rm -rf "+output_folder+"/tmp_netmhc_wt")


def tcr_specificity(prefix, rna_fastq_1, rna_fastq_2, neoantigen_input, mixcr_path, ergo_path, tcr_specificity_software, output_folder, log_folder):
    start_t = time.time()
    mixcr_output_path = os.path.join(output_folder, prefix+"_mixcr_output")
    mixcr_output_file = os.path.join(mixcr_output_path, prefix)
    if not os.path.exists(mixcr_output_path):
        os.mkdir(mixcr_output_path)
    cmd_mixcr = "java -jar "+mixcr_path+" analyze shotgun -s hs --starting-material rna --only-productive \
        --receptor-type tcr \
        "+rna_fastq_1+"  \
        "+rna_fastq_2+" \
        " + mixcr_output_file
    logging.info(cmd_mixcr)
    subprocess.call(cmd_mixcr, shell=True, executable="/bin/bash")
    end_t1 = time.time()
    output_time("Time Takes for TCR repertoire is "+str(round(end_t1-start_t,2))+" s",output_folder)

    cmd_prepare_input = "python rank_software_input.py -m "+mixcr_output_file+" -n "+neoantigen_input+" -o "+output_folder+" -t "+ tcr_specificity_software+" -p "+prefix
    logging.info(cmd_prepare_input)
    subprocess.call(cmd_prepare_input, shell=True, executable="/bin/bash")

    if tcr_specificity_software == "ERGO":
        cmd_ergo = "python "+ergo_path+" mcpas "+output_folder+"/"+prefix+"_cdr_ergo.csv "+output_folder+"/"+prefix+"_tcr_specificity_score.csv"
        logging.info(cmd_ergo)
        subprocess.call(cmd_ergo, shell=True, executable="/bin/bash")     
    else:
        logging.error("There is no indicated neoantigen rank software!")
    cmd_final_output = "python parse_rank_software.py -i "+output_folder+"/"+prefix+"_tcr_specificity_score.csv -n "+neoantigen_input+" -o "+output_folder+" \
        -t "+ tcr_specificity_software+" \
        -p "+prefix
    logging.info(cmd_final_output)
    subprocess.call(cmd_final_output, shell=True, executable="/bin/bash") 

    cmd_add_detail = "python add_detail_info.py -i "+output_folder+"/"+prefix+"_neoantigen_rank_tcr_specificity.tsv \
        -o "+output_folder+" -p "+prefix
    logging.info(cmd_add_detail)
    subprocess.call(cmd_add_detail, shell=True, executable="/bin/bash") 

    end_t2 = time.time()
    output_time("Time Takes for TCR-pMHC Binding Specificity Prediction is "+str(round(end_t2-end_t1,2))+" s",output_folder)    


def main(args_input = sys.argv[1:]):
    input_type =""
    prefix=""
    output_folder=""
    DNA_normal_1=""
    DNA_normal_2=""
    DNA_tumor_1=""
    DNA_tumor_2=""
    RNA_tumor_1=""
    RNA_tumor_2=""
    prioritization_strategy=""
    interval_list=""
    mutation_vcf=""
    HLA_string=""
    """
    Set default parameter thresholds
    """
    binding_affinity = 500
    binding_stability = 1
    tumor_abundance = 1
    agretopicity = 0.1
    foreignness = 1e-12
    tumor_depth = 5
    tumor_vaf = 0.1
    normal_vaf = 0.05
    min_FFPM = 0.1
    USAGE='''
        Main pipeline for neoantigen prediction
        usage: python NeoHunter.py -p <Prefix> -o <OutputFolder> -dn1 <DNANormal1> -dn2 <DNANormal2> \
            -dt1 <DNATumor1> -dt2 <DNATumor2> -rt1 <RNATumor1> -rt2 <RNATumor2> \
            -ps <PrioritizationStrategy> -l <IntervalList> -vcf <MutationVCF> \
            -ba <BindingAffinity> -bs <BindingStability> -ta <TumorAbundance> \
            -ag <Agretopicity> -fo <Foreignness> -td <TumorDepth> -tv <TumorVAF> -na <NormalVAF> -min_FFPM <minFFPM> \
            -hla <HLAString> -t <Thread> at <AlterationType>
            required arguments:
                -p   | --Prefix : Output files prefix
                -o   | --OutputFolder : Output directory
                -dn1 | --DNANormal1 : DNA normal sequence read1 (required if no VCF file as input)
                -dn2 | --DNANormal2 : DNA normal sequence read2 (optional)
                -dt1 | --DNATumor1 : DNA tumor sequence read1 (required if no VCF file as input)
                -dt2 | --DNATumor2 : DNA tumor sequence read2 (optional)
                -rt1 | --RNATumor1 : RNA tumor sequence read1
                -rt2 | --RNATumor2 : RNA tumor sequence read2 (optional)
                -ps  | --PrioritizationStrategy : Type of prioritization strategy
            optional arguments:
                -l   | --IntervalList: bed file for subsets of genomics regions (optional)
                -vcf | --MutationVCF : VCF file of detected mutations (optional)
                -ba  | --BindingAffinity : peptide-MHC binding affinity threshold (optional, default is 500)
                -bs  | --BindingStability : peptide-MHC binding stability threshold (optional, default is 1.4)
                -ta  | --TumorAbundance : tumor abundance threshold (optional, default is 6)
                -ag  | --Agretopicity : Agretopicity threshold (optional, default is 0.1)
                -fo  | --Foreignness : Foreignness threshold (optional, default is 10e-16)
                -td  | --TumorDepth : tumor allele depth threshold (optional, default is 5)
                -tv  | --TumorVAF : tumor variant allele frequency threshold (optional, default is 0.1)
                -nv  | --NormalVAF : normal variant allele frequency threshold (optional, default is 0.05)
                -min_FFPM |--minFFPM : FFPM thredhold for STAR-Fusion (optional, default is 0.1)
                -hla | --HLAString : Four-digit HLA string input, in the format like HLA-A01:01 (optional)
                -t | --Thread : Number of thread running in parallel (optional, default is 8)
                -at  | --AlterationType : type of mutation to find neoantigen (optional, default is : snv, indel, fusion, splicing)
    '''
    
    parser = argparse.ArgumentParser()
    # parser.add_argument(
    #     '-i','--InputType',
    #     help="Type of input files (DNA/VCF)"
    # )
    parser.add_argument(
        '-p','--Prefix',
        help="Output files prefix"
    )
    parser.add_argument(
        '-o','--OutputFolder',
        help="Output directory"
    )
    parser.add_argument(
        '-dn1','--DNANormal1',
        help="DNA normal sequence read1"
    )
    parser.add_argument(
        '-dn2','--DNANormal2',
        help="DNA normal sequence read2 (optional)"
    )
    parser.add_argument(
        '-dt1','--DNATumor1',
        help="DNA tumor sequence read1"
    )
    parser.add_argument(
        '-dt2','--DNATumor2',
        help="DNA tumor sequence read2 (optional)"
    )
    parser.add_argument(
        '-rt1','--RNATumor1',
        help="RNA tumor sequence read1"
    )
    parser.add_argument(
        '-rt2','--RNATumor2',
        help="RNA tumor sequence read2 (optional)"
    )
    parser.add_argument(
        '-ps','--PrioritizationStrategy',
        help="Type of prioritization strategy"
    )
    parser.add_argument(
        '-l','--IntervalList',
        help="subsets of genomic regions (optional)"
    )
    parser.add_argument(
        '-vcf','--MutationVCF',
        help="VCF file of detected mutations"
    )
    parser.add_argument(
        '-ba','--BindingAffinity',
        help="peptide-MHC binding affinity threshold (optional, default is 500)"
    )
    parser.add_argument(
        '-bs','--BindingStability',
        help="peptide-MHC binding stability threshold (optional, default is 1.4)"
    )
    parser.add_argument(
        '-ta','--TumorAbundance',
        help="tumor abundance threshold (optional, default is 6)"
    )
    parser.add_argument(
        '-ag','--Agretopicity',
        help="Agretopicity threshold (optional, default is 0.1)"
    )
    parser.add_argument(
        '-fo','--Foreignness',
        help="Foreignness threshold (optional, default is 10e-16)"
    )
    parser.add_argument(
        '-td','--TumorDepth',
        help="tumor allele depth threshold (optional, default is 5)"
    )
    parser.add_argument(
        '-tv','--TumorVAF',
        help="tumor variant allele frequency threshold (optional, default is 0.1)"
    )
    parser.add_argument(
        '-nv','--NormalVAF',
        help="normal variant allele frequency threshold (optional, default is 0.05)"
    )
    parser.add_argument(
        '-min_FFPM','--minFFPM',
        help="FFPM thredhold for STAR-Fusion (optional, default if 0.1)"
    )
    parser.add_argument(
        '-hla','--HLAString',
        help="Four-digit HLA string input, in the format like HLA-A01:01 (optional)"
    )
    parser.add_argument(
        '-t','--Thread',
        help="Number of thread running in parallel (optional)"
    )
    parser.add_argument(
        '-at','--AlterationType',
        help="type of mutation to find neoantigen, default is all: SNV, INDEL, FUSION, SPLICING (optional)"
    )
    args = parser.parse_args(args_input)
    # input_type=args.InputType
    if not args.OutputFolder:
        logging.error("Output folder is required")
        sys.exit(2)
    if not args.RNATumor1:
        logging.error("RNA tumor read1 is required")
        sys.exit(2)
    output_folder=args.OutputFolder
    RNA_tumor_1=args.RNATumor1
    HLA_string=args.HLAString
    prioritization_strategy=args.PrioritizationStrategy
    thread=8
    alteration_type="snv,indel,fusion,splicing"
    if (args.Prefix):
        if (args.Prefix != "None"):
            prefix = args.Prefix
    if (args.DNANormal2):
        if (args.DNANormal2 != "None"):
            DNA_normal_2=args.DNANormal2
    if (args.DNATumor2):
        if (args.DNATumor2 != "None"):
            DNA_tumor_2=args.DNATumor2
    if (args.RNATumor2):
        if (args.RNATumor2 != "None"):
            RNA_tumor_2=args.RNATumor2
    if (args.MutationVCF):
        if (args.MutationVCF != "None"):
            mutation_vcf=args.MutationVCF
    # Get parameter thresholds
    if (args.BindingAffinity):
        if (args.BindingAffinity != "None"):
            binding_affinity = float(args.BindingAffinity)
    if (args.BindingStability):
        if (args.BindingStability != "None"):
            binding_stability=float(args.BindingStability)
    if (args.TumorAbundance):
        if (args.TumorAbundance != "None"):
            tumor_abundance=float(args.TumorAbundance)
    if (args.Agretopicity):
        if (args.Agretopicity != "None"):
            agretopicity=float(args.Agretopicity)
    if (args.Foreignness):
        if (args.Foreignness != "None"):
            foreignness=float(args.Foreignness)
    if (args.TumorDepth):
        if (args.TumorDepth != "None"):
            tumor_depth=int(args.TumorDepth)
    if (args.TumorVAF):
        if (args.TumorVAF != "None"):
            tumor_vaf=float(args.TumorVAF)
    if (args.NormalVAF):
        if (args.NormalVAF != "None"):
            normal_vaf=float(args.NormalVAF)
    if (args.minFFPM):
        if (args.minFFPM != "None"):
            min_FFPM=float(args.minFFPM)
    if (args.Thread):
        if (args.Thread != "None"):
            thread = int(args.Thread)
    if (args.IntervalList):
        if (args.IntervalList != "None"):
            interval_list=args.IntervalList
    if (args.AlterationType):
        if (args.AlterationType != "None"):
            alteration_type=args.AlterationType
    if (output_folder ==""):
        logging.info(USAGE)
        sys.exit(2)
    if (mutation_vcf == ""):
        DNA_normal_1=args.DNANormal1
        DNA_tumor_1=args.DNATumor1
        if (DNA_normal_1 =="" or DNA_tumor_1 =="" or DNA_normal_1 =="None" or DNA_tumor_1 =="None"):
            logging.error("Please provide all necessary DNA files")
            logging.info(USAGE)
            sys.exit(2)
        if (RNA_tumor_1 =="" or RNA_tumor_1 =="None"):
            logging.error("Please provide at least one RNA files")
            logging.info(USAGE)
            sys.exit(2)
    supported_prioritization_strategies=["direct", "indirect", "both"]
    if prioritization_strategy not in supported_prioritization_strategies:
        logging.error("Unrecognized Prioritization Strategy. Please provide one of the following: " + str(supported_prioritization_strategies))
        sys.exit(2)
    
    if os.path.exists(output_folder):
        logging.info("Output folder already exists")
    else:
        os.mkdir(output_folder)

    config_file = "config.yaml"
    f=open(config_file)
    config_list=yaml.safe_load(f)
    """
    Parse dataset in config file
    """
    mills_path="/data8t_2/zzt/data/variant_call/pTuneos/database/VCF_annotation/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" 
    # "./database/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
    dbsnp_path="/data8t_5/zzt/reference_annotation/dbsnp_138.hg19.vcf"
    # "./database/dbsnp_138.hg19.vcf.gz"
    OneKG_path="/data8t_5/zzt/reference_annotation/1000G_phase1.indels.hg19.sites.vcf"
    # "./database/1000G_phase1.indels.hg19.sites.vcf.gz"
    human_peptide_path="/data8t_2/zzt/data/variant_call/neo_pipeline/database/Homo_sapiens.GRCh37.pep.all.fa"
    # "./database/Homo_sapiens.GRCh37.pep.all.fa"
    kallisto_cdna_path="/data8t_2/zzt/data/variant_call/neo_pipeline/database/Homo_sapiens.GRCh37.cdna.all.fa"
    # "./database/Homo_sapiens.GRCh37.cdna.all.fa"
    funcotator_source_path="/data8t_2/zzt/data/variant_call/pTuneos/database/funcotator_dataSources.v1.6.20190124g"
    # "./database/funcotator_dataSources.v1.6.20190124g"
    star_fusion_dataset_path = "/data8t_5/zzt/reference_annotation/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/"
    # "./database/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/"
    iedb_path="/data8t_2/zzt/antigen.garnish/iedb.fasta"
    # "./database/iedb.fasta"
    reference="/data8t_2/zzt/data/variant_call/neo_pipeline/database/hg19.fa" 
    # "./database/hg19.fa"
    refseq_ann="/data8t_2/zzt/data/variant_call/neo_pipeline/database/hg19.refGene.gtf.gz"
    # "./database/hg19.refGene.gtf.gz"
    if config_list["mills_path"] != "None":
        mills_path = config_list["mills_path"]
    if config_list["dbsnp_path"] != "None":
        dbsnp_path = config_list["dbsnp_path"]
    if config_list["OneKG_path"] != "None":
        OneKG_path = config_list["OneKG_path"]
    if config_list["human_peptide_path"] != "None":
        human_peptide_path = config_list["human_peptide_path"]
    if config_list["kallisto_cdna_path"] != "None":
        kallisto_cdna_path = config_list["kallisto_cdna_path"]
    if config_list["funcotator_source_path"] != "None":
        funcotator_source_path = config_list["funcotator_source_path"]
    if config_list["star_fusion_dataset_path"] != "None":
        star_fusion_dataset_path = config_list["star_fusion_dataset_path"]
    if config_list["iedb_path"] != "None":
        iedb_path = config_list["iedb_path"]
    if config_list["reference"] != "None":
        reference = config_list["reference"]
    if config_list["refseq_ann"] != "None":
        refseq_ann = config_list["refseq_ann"]
    """
    Parse software path in config file
    """
    annotation_software = config_list["annotation_software"] # default is VEP
    tcr_specificity_software = config_list["tcr_specificity_software"] # default if ERGO
    asneo_path = config_list["asneo_path"]
    bwa_path = config_list["bwa_path"]
    ergo_path = config_list["ergo_path"]
    gatk_path = config_list["gatk_path"]
    kallisto_path = config_list["kallisto_path"]
    netmhc_path = config_list["netmhc_path"]
    netMHCstabpan_path = config_list["netMHCstabpan_path"]
    mixcr_path = config_list["mixcr_path"]
    optitype_path = config_list["optitype_path"]
    perl_path = config_list["perl_path"]
    picard_path = config_list["picard_path"]
    samtools_path = config_list["samtools_path"]
    snpeff_path = config_list["snpeff_path"]
    star_fusion_path = config_list["star_fusion_path"]
    star_path = config_list["star_path"]
    vep_path = config_list["vep_path"]
    vep_cache = config_list["vep_cache"]

    """
    Validate dataset and software paths
    """
    if not os.path.exists(bwa_path):
        logging.error("BWA does not exist")
        os._exit(1)
    if not os.path.exists(gatk_path):
        logging.error("GATK does not exist")
        os._exit(1)
    if not os.path.exists(kallisto_path):
        logging.error("Kallisto does not exist")
        os._exit(1)
    if not os.path.exists(netmhc_path):
        logging.error("Netmhc does not exist")
        os._exit(1)
    if not os.path.exists(optitype_path):
        logging.error("Optitype does not exist")
        os._exit(1)
    if not os.path.exists(picard_path):
        logging.error("Picard.jar does not exist")
        os._exit(1)
    if not os.path.exists(samtools_path):
        logging.error("Samtools does not exist")
        os._exit(1)
    if not os.path.exists(snpeff_path):
        logging.error("SnpEff does not exist")
        os._exit(1)
    if not os.path.exists(vep_path):
        logging.error("VEP does not exist")
        os._exit(1)
    if not os.path.exists(vep_cache + "/homo_sapiens"):
        logging.error("VEP cache (homo_sapiens) does not exist")
        os._exit(1)

    log_folder = os.path.join(output_folder,"log")
    if not os.path.exists(log_folder):
        os.mkdir(log_folder)

    info_folder = os.path.join(output_folder,"info")
    if not os.path.exists(info_folder):
        os.mkdir(info_folder)

    time1 = time.time()
    cpu_time_start = psutil.cpu_times()
    logging.info("[PART1] Start Alignment & Alteration Detection & HLA Typing")
    if mutation_vcf != "":
        logging.info("Skip Alignment and SNV/Indel Detection")
    else:
        logging.info("Start Alignment") 
        alignment_folder = os.path.join(output_folder,"genome_alignment")
        if not os.path.exists(alignment_folder):
            os.mkdir(alignment_folder)
        alteration_detection_folder = os.path.join(output_folder,"alteration_detection")
        if not os.path.exists(alteration_detection_folder):
            os.mkdir(alteration_detection_folder)
        DNA_normal_bam = ""
        DNA_tumor_bam = ""
        convert_process = []
        if (DNA_normal_1.endswith("fastq") or DNA_normal_1.endswith("fastq.gz")):
            convert_normal=multiprocessing.Process(target=fastq2bam,args=(prefix+"_normal", DNA_normal_1, DNA_normal_2, alignment_folder, reference, bwa_path, thread))
            convert_process.append(convert_normal)
            DNA_normal_bam = alignment_folder+"/tmp_"+prefix+"_normal.sam"
        else:
            DNA_normal_bam = DNA_normal_1
        if (DNA_tumor_1.endswith("fastq") or DNA_tumor_1.endswith("fastq.gz")):
            convert_tumor=multiprocessing.Process(target=fastq2bam,args=(prefix+"_tumor", DNA_tumor_1, DNA_tumor_2, alignment_folder, reference, bwa_path, thread))
            convert_process.append(convert_tumor)
            DNA_tumor_bam = alignment_folder+"/tmp_"+prefix+"_tumor.sam"
        else:
            DNA_tumor_bam = DNA_tumor_1
        for p in convert_process:
            p.daemon = True
            p.start()
        for p in convert_process:
            p.join()
        logging.info("Finish Alignment")
        
        alignment_process = []
        normal_process=multiprocessing.Process(target=alignment,args=(prefix+"_normal", DNA_normal_bam, alignment_folder, log_folder, samtools_path, thread))
        alignment_process.append(normal_process)
        tumor_process=multiprocessing.Process(target=alignment,args=(prefix+"_tumor", DNA_tumor_bam, alignment_folder, log_folder, samtools_path, thread))
        alignment_process.append(tumor_process)
        for p in alignment_process:
            p.daemon = True
            p.start()
        for p in alignment_process:
            p.join()

        snv_indel_folder = os.path.join(alteration_detection_folder,prefix+"_snv_indel_output")
        if not os.path.exists(snv_indel_folder):
            os.mkdir(snv_indel_folder)

        pre_mutect2_process = []
        normal_process=multiprocessing.Process(target=pre_mutect2,args=(prefix+"_normal", reference, alignment_folder, snv_indel_folder, log_folder, OneKG_path, mills_path, dbsnp_path, gatk_path, picard_path, samtools_path))
        pre_mutect2_process.append(normal_process)
        tumor_process=multiprocessing.Process(target=pre_mutect2,args=(prefix+"_tumor", reference, alignment_folder, snv_indel_folder, log_folder, OneKG_path, mills_path, dbsnp_path, gatk_path, picard_path, samtools_path))
        pre_mutect2_process.append(tumor_process)
        for p in pre_mutect2_process:
            p.daemon = True
            p.start()
        for p in pre_mutect2_process:
            p.join()
        time2 = time.time()
        output_time("Time Takes for Preparing SNV/Indel Detection is "+str(round(time2-time1,2))+" s",alteration_detection_folder)

        logging.info("Detect SNV/Indel with GATK")
        mutation_detection(prefix, reference, snv_indel_folder, log_folder, gatk_path, interval_list, tumor_depth, tumor_vaf, normal_vaf)
        logging.info("Detect SNV/Indel with GATK finished")
        time3 = time.time()
        output_time("Time Takes for SNV/Indel Detection is "+str(round(time3-time2,2))+" s",alteration_detection_folder)

    rna_fastq_path_1 = ""
    rna_fastq_path_2 = ""
    if (RNA_tumor_1.endswith("bam")):
        cmd_convert_rna_1 = "java -Xmx4G -jar " + picard_path + " SamToFastq I="+RNA_tumor_1+" FASTQ="+output_folder+"/RNA_tumor_1.fastq"
        logging.info(cmd_convert_rna_1)
        os.system(cmd_convert_rna_1)
        cmd_convert_rna_2 = "java -Xmx4G -jar " + picard_path + " SamToFastq I="+RNA_tumor_2+" FASTQ="+output_folder+"/RNA_tumor_2.fastq"
        logging.info(cmd_convert_rna_2)
        os.system(cmd_convert_rna_2)
        rna_fastq_path_1 = output_folder+"/RNA_tumor_1.fastq"
        rna_fastq_path_2 = output_folder+"/RNA_tumor_2.fastq"
    else:
        rna_fastq_path_1 = RNA_tumor_1
        rna_fastq_path_2 = RNA_tumor_2

    time4 = time.time()

    logging.info("Start HLA Typing and RNA-seq Analysis")
    prepare_data_process = []

    hla_type = multiprocessing.Process(target=hla_typing,args=(prefix, rna_fastq_path_1, rna_fastq_path_2, output_folder, log_folder, optitype_path))
    prepare_data_process.append(hla_type)

    transcript_quantification_folder = os.path.join(alteration_detection_folder,prefix+"_transcript_quantification")
    if not os.path.exists(transcript_quantification_folder):
        os.mkdir(transcript_quantification_folder)
    rna_process = multiprocessing.Process(target=kallisto_expression,args=(rna_fastq_path_1,rna_fastq_path_2,kallisto_cdna_path,kallisto_path,transcript_quantification_folder,log_folder,thread))
    prepare_data_process.append(rna_process)
    for p in prepare_data_process:
        p.daemon = True
        p.start()
    for p in prepare_data_process:
        p.join()
    logging.info("Finish HLA Typing and RNA-seq Analysis")

    time5 = time.time()
    output_time("Time Takes for data preprocessing is "+str(round(time5-time4,2))+" s",alteration_detection_folder)

    hla_str=""
    if HLA_string != "None" and HLA_string != None:
        hla_str = HLA_string
    else:
        # read optitype output to string
        hla_output_path = os.path.join(output_folder,"hla_typing")
        hla_typing_file=hla_output_path+"/"+prefix+"_hla.tsv"
        with open(hla_typing_file) as file:
            reader = csv.reader(file, delimiter="\t")
            output_hla_string = ""
            for line in reader:
                if line[0] == "":
                    continue
                for i in range(1,7,1):
                    output_hla_string += "HLA-"
                    output_hla_string += line[i].replace("*","")
                    output_hla_string += ","
            hla_str = output_hla_string[:-1]
    logging.info("[PART2] Start Peptide-MHC Binding Prediction")
    pmhc_binding_prediction_folder = os.path.join(output_folder,"pMHC_binding_prediction")
    if not os.path.exists(pmhc_binding_prediction_folder):
        os.mkdir(pmhc_binding_prediction_folder)
    if (mutation_vcf == ""or mutation_vcf=="None"):
        mutation_vcf = os.path.join(snv_indel_folder, prefix+"_filter.vcf")#output_folder+"/filter.vcf"

    snv_indel_pred(prefix, mutation_vcf, reference, human_peptide_path, tumor_abundance, alteration_detection_folder, perl_path, gatk_path, vep_path, vep_cache, snpeff_path, funcotator_source_path, annotation_software)

    if ("fusion" in alteration_type) or ("splicing" in alteration_type):
        fusion_pred(prefix, RNA_tumor_1, RNA_tumor_2, star_fusion_dataset_path, alteration_detection_folder, star_fusion_path, star_path, tumor_abundance, min_FFPM, thread)
    else:
        logging.info("[WARNING] Skipping fusion mutation neoantigen prediction")
    if ("splicing" in alteration_type):
        splicing_pred(prefix, RNA_tumor_1,RNA_tumor_2,refseq_ann,star_path,star_fusion_dataset_path,asneo_path,alteration_detection_folder,reference, tumor_abundance, thread)
    else:
        logging.info("[WARNING] Skipping splicing mutation neoantigen prediction")
    os.system("cp "+snv_indel_folder+"/"+prefix+"_snv_indel_annotation.vcf "+info_folder+"/"+prefix+"_snv_indel.vcf")
    os.system("cp "+alteration_detection_folder+"/"+prefix+"_fusion_output/star-fusion.fusion_predictions.abridged.coding_effect.tsv "+info_folder+"/"+prefix+"_fusion.tsv")
    logging.info("Finish Mutation Annotation, and Convert to Fasta")

    logging.info("Concat fasta files from different mutation types")
    tmp_fasta_folder = os.path.join(alteration_detection_folder,prefix+"_tmp_fasta")
    cmd_cat = "cat "+alteration_detection_folder+"/"+prefix+"_snv_indel.fasta "+alteration_detection_folder+"/"+prefix+"_fusion.fasta > "+tmp_fasta_folder+"/"+prefix+"_snv_indel_fusion.fasta"
    logging.info(cmd_cat)
    os.system(cmd_cat)

    time6 = time.time()
    logging.info("Start netmhcpan peptide-mhc binding affinity prediction")
    logging.info(hla_str)
    mutation_netmhc_parallel(prefix, netmhc_path, tmp_fasta_folder, pmhc_binding_prediction_folder, hla_str)
    logging.info("Finish netmhcpan peptide-mhc binding affinity prediction")
    
    time7 = time.time()
    output_time("Time Takes for Binding Affinity Prediction is "+str(round(time7-time6,2))+" s",pmhc_binding_prediction_folder)

    cmd_cat = "cat "+tmp_fasta_folder+"/"+prefix+"_snv_indel_fusion.fasta "+alteration_detection_folder+"/"+prefix+"_splicing.fasta > "+alteration_detection_folder+"/"+prefix+"_alteration_derived_pep.fasta"
    logging.info(cmd_cat)
    os.system(cmd_cat)
    parse_netmhc_snv_indel = "python parse_netMHC.py -i "+pmhc_binding_prediction_folder+" \
        -g "+alteration_detection_folder+"/"+prefix+"_alteration_derived_pep.fasta -o "+pmhc_binding_prediction_folder+" \
        -b "+str(binding_affinity)+" -l " + hla_str+" -p "+prefix
    logging.info(parse_netmhc_snv_indel)
    subprocess.call(parse_netmhc_snv_indel, shell=True, executable="/bin/bash")

    time8 = time.time()
    logging.info("Start netmhcstabpan peptide-mhc binding stability prediction")
    run_calculation = "python bindstab_filter.py -i " + pmhc_binding_prediction_folder+"/"+prefix+"_bindaff_filtered.tsv \
            -o "+pmhc_binding_prediction_folder+" -n "+netMHCstabpan_path+" \
            -b "+str(binding_stability) + " -p "+prefix
    logging.info(run_calculation)
    subprocess.call(run_calculation, shell=True, executable="/bin/bash")
    logging.info("Finish netmhcstabpan peptide-mhc binding stability prediction")
    time9 = time.time()
    output_time("Time Takes for Binding Stability Prediction is "+str(round(time9-time8,2))+" s",pmhc_binding_prediction_folder)
 
    logging.info("[PART3] Start Immunogenicity Evaluation & Neoantigen Prioritization")
    prioritization_folder = os.path.join(output_folder,"prioritization")
    alteration_type = alteration_type.replace(" ", "")
    if not os.path.exists(prioritization_folder):
        os.mkdir(prioritization_folder)
    if prioritization_strategy=="indirect" or prioritization_strategy=="both" :
        time10 = time.time()
        cmd_recognition_associated_prioritization = "python recognition_associated_prioritization.py -i " + pmhc_binding_prediction_folder+" -I "+iedb_path+"\
        -o "+prioritization_folder+" -n "+netmhc_path+" -a "+str(agretopicity)+" -f "+str(foreignness)+" -t "+alteration_type+" -p "+prefix
        logging.info(cmd_recognition_associated_prioritization)
        subprocess.call(cmd_recognition_associated_prioritization, shell=True, executable="/bin/bash")
        time11 = time.time()
        output_time("Time Takes for Agretopicity&Foreignness Feature Evaluation & Binding-affinity-related Prioritization is "+str(round(time11-time10,2))+" s",prioritization_folder)
    if prioritization_strategy=="direct" or prioritization_strategy=="both":
        tcr_specificity(prefix, rna_fastq_path_1, rna_fastq_path_2, pmhc_binding_prediction_folder+"/"+prefix+"_candidate_pmhc.csv", mixcr_path, ergo_path, tcr_specificity_software, prioritization_folder, log_folder)
    
    time10 = time.time()
    cpu_time_end = psutil.cpu_times()
    start_time = float(str(cpu_time_start).strip().split(',')[0].split('=')[1])
    end_time = float(str(cpu_time_end).strip().split(',')[0].split('=')[1])
    output_time("Time Takes for NeoHunter Pipeline is "+str(round(time10-time1,2))+" s",alteration_detection_folder)
    output_time("Total CPU Time for NeoHunter Pipeline is "+str(round(end_time-start_time,2))+" s",alteration_detection_folder)    
    logging.info("Finish Immunogenicity Evaluation & Neoantigen Prioritization, please check neoantigen_rank.csv in output folder")


if __name__ == '__main__':
    main()