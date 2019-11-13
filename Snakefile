import os
import paramiko
import glob
import re
from collections import Counter
import fnmatch
import pandas as pd

'''
The file hierarchy looks like:
    
    ./HorseGeneAnnotation/
        private/
            sequence/
                RNASEQ/
		    bam/
		        {NCBI_GCF,ENSEMBL_GCA}/
		            {single_end,paired_end}/
        public/
	    refgen/
	        {NCBI_GCF,ENSEMBL_GCA}/
		    STAR_INDICES/
       	            	{single_end,paired_end}/
	                    GFF/
		                Merged/

'''

from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
s3_key_id = os.environ.get('AWS_ACCESS_KEY')
s3_access_key = os.environ.get('AWS_SECRET_KEY')

S3 = S3RemoteProvider(
    endpoint_url='https://s3.msi.umn.edu', 
    access_key_id=s3_key_id, 
    secret_access_key=s3_access_key
)

configfile: "config.yaml" 

from snakemake.remote.SFTP import RemoteProvider as SFTPRemoteProvider
key = paramiko.agent.Agent().get_keys()[0]

#SFTP = SFTPRemoteProvider('login.msi.umn.edu',username='cull0084',private_key='/project/cull0084/.ssh/id_rsa')
SFTP = SFTPRemoteProvider(username='cull0084',private_key=key)

#se_samples_tmp = []
#for path in [x for x in S3._s3c.list_keys('HorseGeneAnnotation') if 'fastq' in x]: 
#    dir_path, se_file = os.path.split(path)
#    se_file = re.search(r'^(.*?)(_R[1-2]_)',se_file).group(1)
#    se_samples_tmp.append(se_file)
#    
#SE_SAMPLES = [k for k, v in Counter(se_samples_tmp).items() if v == 1]

#SUBSET = ["61M", "35M", "87M", "52M", "17M", "36F", "90M", "66M", "67M", "37M"]
SUBSET = ['61M', '35M', '87M', '52M', '17M', '36F', '90M', '66M', '67M', '37M',
          '65F', '67F', '82M', '69F', '12F', '20M', '86M', '39F', '1M', '49M', 
          '66F', '65M', '35F', '51F', '40M', '39M']



#hiseq
#ALL_SAMPLES, = SFTP.glob_wildcards("login.msi.umn.edu/panfs/roc/data_release/3/umgc/pre2018/2015-q4/mccuem/hiseq/151006_D00635_0082_BC7HAHANXX/Project_McCue_Project_022/{sample}_R2_001.fastq")

#novaseq
ALL_SAMPLES, = SFTP.glob_wildcards("login.msi.umn.edu//panfs/roc/data_release/3/umgc/mccuem/novaseq/181112_A00223_0050_AHCWF7DSXX/McCue_Project_032/{sample}_R2_001.fastq.gz")


SAMPLES = []
for i in SUBSET:
    i = i[:-1] + "_" + i[-1]
    for j in ALL_SAMPLES:
        if i in j:
            SAMPLES.append(j)

print(SAMPLES)

#REF_GFF = [f"{config['NCBI_GCF']}", f"{config['ENSEMBL_GCA']}"]
REF_GFF = [f"{config['NCBI_GCF']}"]
#REF_GFF = [f"{config['ENSEMBL_GCA']}"]

rule all:
    input:
#        expand('HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/Merged/{sample}/{sample}.gff',sample=SAMPLES,GCF=REF_GFF),
#        f"HorseGeneAnnotation/public/refgen/{config['NCBI_GCF']}/paired_end/GFF/Merged/transcript_fpkm.tsv",
#        f"HorseGeneAnnotation/public/refgen/{config['NCBI_GCF']}/paired_end/GFF/Merged/gene_fpkm.tsv"
##        expand('HorseGeneAnnotation/private/sequence/RNASEQ/bam/{GCF}/paired_end/{sample}_Aligned.out.bam',GCF=REF_GFF,sample=SAMPLES)
##        expand('HorseGeneAnnotation/private/sequence/RNASEQ/bam/{GCF}/paired_end/{sample}.sorted.bam',GCF=REF_GFF,sample=SAMPLES)
##        expand('HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/{sample}.gff',GCF=REF_GFF,sample=SAMPLES)
##        expand('HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/Merged.gff',GCF=REF_GFF),
##       expand('HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/Merged/{sample}/{sample}.gff',sample=SAMPLES,GCF=REF_GFF)
        expand('HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/Merged/transcript_fpkm.tsv',GCF=REF_GFF),
        expand('HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/Merged/gene_fpkm.tsv',GCF=REF_GFF),
        expand('HorseGeneAnnotation/public/refgen/{GCF}/paired_end/SalmonQuant/novaseq_TPM.tsv',GCF=REF_GFF)

# ----------------------------------------------------------
#       Trimming
# ----------------------------------------------------------

rule pe_trim_reads:
    input:
        #R1 = S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}_R1_001.fastq.gz'),
        #R2 = S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}_R2_001.fastq.gz')
        R1 = SFTP.remote('login.msi.umn.edu/panfs/roc/data_release/3/umgc/mccuem/novaseq/181112_A00223_0050_AHCWF7DSXX/McCue_Project_032/{sample}_R1_001.fastq.gz'),
        R2 = SFTP.remote('login.msi.umn.edu/panfs/roc/data_release/3/umgc/mccuem/novaseq/181112_A00223_0050_AHCWF7DSXX/McCue_Project_032/{sample}_R2_001.fastq.gz')
    output:
        R1 = 'HorseGeneAnnotation/private/sequence/trimmed/{sample}_trim1.fastq.gz',
        R2 = 'HorseGeneAnnotation/private/sequence/trimmed/{sample}_trim2.fastq.gz'

    message:
        'AdapterRemoval - removing adapters and low quality bases on {wildcards.sample}'
    shell:
        '''
        AdapterRemoval \
        --file1 {input.R1} \
        --file2 {input.R2} \
        --output1 {output.R1} \
        --output2 {output.R2} \
        --gzip \
        --trimns \
        --trimqualities \
        '''

#rule get_se_fastqs:
#    input:
#        S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}_R1_001.fastq.gz')
#    output:
#        'local/HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}_R1_001.fastq.gz'
#    run:
#        shell('cp {input[0]} {output[0]}')


#rule se_trim_read:
#    input:
#        #R1 = '/scratch/single_end_mapping/se_fastq/{sample}_R1_001.fastq.gz' 
#        R1 = S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}_R1_001.fastq.gz',keep_local=True)
#    output:
#        #R1 = '/scratch/single_end_mapping/trimmed_data/{sample}_se_trim.fastq.gz'
#        R1 = temp('se_trimmed_data/{sample}_trim1.fastq.gz')
#    message:
#        'AdapterRemoval - removing adapters and low quality bases on SE reads {wildcards.sample}'
#    shell:
#        '''
#        AdapterRemoval \
#        --file1 {input.R1} \
#        --output1 {output.R1} \
#        --gzip \
#        --trimns \
#        --trimqualities \
#        --minquality 10 \
#        '''

# ----------------------------------------------------------
#       QC
# ----------------------------------------------------------

#rule qc_trim:
#    input:
#        S3.remote('trimmed_data/{sample}.fastq.gz')
#    params:
#        out_dir = S3.remote('qc/qc_trim')
#    output:
#        S3.remote('qc/qc_trim/{sample}_fastqc.html')
#    message:
#        'FastQC - performing quality control on trimmed {wildcards.sample}'
#    shell:
#        '''
#        fastqc \
#        -o {params.out_dir} \
#        -f fastq \
#        {input}
#        '''
#
#rule qc_raw:
#    input:
#        S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}.fastq.gz')
#    params:
#        out_dir = S3.remote('qc/qc_raw')
#    output:
#        S3.remote('qc/qc_raw/{sample}_fastqc.html')
#    message:
#        'FastQC - performing quality control on {wildcards.sample}'
#    shell:
#        '''
#        fastqc \
#        -o {params.out_dir} \
#        -f fastq \
#        {input}
#        '''

# ----------------------------------------------------------
#       STAR indices
# ----------------------------------------------------------

rule download_STAR:
    input:
        expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/download.done',GCF=config['NCBI_GCF']),
        expand('HorseGeneAnnotation/public/refgen/{GCA}/STAR_INDICES/download.done',GCA=config['ENSEMBL_GCA'])

#DOES NOT STORE UNMAPPED READS (--outReadsUnmapped) nor log output
rule ncbi_STAR_index:
    input:
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/Genome',GCF=config['NCBI_GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/SA',GCF=config['NCBI_GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/SAindex',GCF=config['NCBI_GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/chrLength.txt',GCF=config['NCBI_GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/chrName.txt',GCF=config['NCBI_GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/chrNameLength.txt',GCF=config['NCBI_GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/chrStart.txt',GCF=config['NCBI_GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/exonGeTrInfo.tab',GCF=config['NCBI_GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/exonInfo.tab',GCF=config['NCBI_GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/geneInfo.tab',GCF=config['NCBI_GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/genomeParameters.txt',GCF=config['NCBI_GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/sjdbInfo.txt',GCF=config['NCBI_GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/sjdbList.fromGTF.out.tab',GCF=config['NCBI_GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/sjdbList.out.tab',GCF=config['NCBI_GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/transcriptInfo.tab',GCF=config['NCBI_GCF']),keep_local=True) 
    output:
        touch(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/download.done',GCF=config['NCBI_GCF']))


rule ensembl_STAR_index:
    input:
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCA}/STAR_INDICES/Genome',GCA=config['ENSEMBL_GCA']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCA}/STAR_INDICES/SA',GCA=config['ENSEMBL_GCA']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCA}/STAR_INDICES/SAindex',GCA=config['ENSEMBL_GCA']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCA}/STAR_INDICES/chrLength.txt',GCA=config['ENSEMBL_GCA']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCA}/STAR_INDICES/chrName.txt',GCA=config['ENSEMBL_GCA']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCA}/STAR_INDICES/chrNameLength.txt',GCA=config['ENSEMBL_GCA']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCA}/STAR_INDICES/chrStart.txt',GCA=config['ENSEMBL_GCA']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCA}/STAR_INDICES/exonGeTrInfo.tab',GCA=config['ENSEMBL_GCA']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCA}/STAR_INDICES/exonInfo.tab',GCA=config['ENSEMBL_GCA']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCA}/STAR_INDICES/geneInfo.tab',GCA=config['ENSEMBL_GCA']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCA}/STAR_INDICES/genomeParameters.txt',GCA=config['ENSEMBL_GCA']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCA}/STAR_INDICES/sjdbInfo.txt',GCA=config['ENSEMBL_GCA']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCA}/STAR_INDICES/sjdbList.fromGTF.out.tab',GCA=config['ENSEMBL_GCA']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCA}/STAR_INDICES/sjdbList.out.tab',GCA=config['ENSEMBL_GCA']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCA}/STAR_INDICES/transcriptInfo.tab',GCA=config['ENSEMBL_GCA']),keep_local=True) 
    output:
        touch(expand('HorseGeneAnnotation/public/refgen/{GCA}/STAR_INDICES/download.done',GCA=config['ENSEMBL_GCA']))

# ----------------------------------------------------------
#       STAR mapping
# ----------------------------------------------------------

rule pe_STAR_mapping:
    input:
        R1 = 'HorseGeneAnnotation/private/sequence/trimmed/{sample}_trim1.fastq.gz',
        R2 = 'HorseGeneAnnotation/private/sequence/trimmed/{sample}_trim2.fastq.gz',
        star_dl = 'HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/download.done'
    output:
        S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/bam/{GCF}/paired_end/{sample}_Aligned.out.bam')
    params:
        out_prefix = S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/bam/{GCF}/paired_end/{sample}_'),
        star_index = 'HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES'
    message:
        'STAR - Creating: {output} '
    run:
        assert os.path.exists(input.star_dl)
        shell('''
        STAR \
        --genomeDir {params.star_index} \
        --genomeLoad LoadAndKeep \
        --readFilesIn {input.R1} {input.R2} \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix {params.out_prefix} \
        --outSAMtype BAM Unsorted \
        ''')

##rule pe_STAR_mapping:
##    input:
##        R1 = 'pe_trimmed_data/{sample}_trim1.fastq.gz',
##        R2 = 'pe_trimmed_data/{sample}_trim2.fastq.gz',
##        star_index = expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/download.done',GCF=config['NCBI_GCF'])
##    output:
##        S3.remote(expand('HorseGeneAnnotation/private/sequence/RNASEQ/bam/{GCF}/paired_end/{sample}_Aligned.out.bam',sample=SAMPLES,GCF=config['NCBI_GCF']),keep_local=True)
##    params:
##        out_prefix = S3.remote(expand('HorseGeneAnnotation/private/sequence/RNASEQ/bam/{GCF}/paired_end/{sample}_',sample=SAMPLES,GCF=config['NCBI_GCF']),keep_local=True),
##        star_index = expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES',GCF=config['NCBI_GCF'])
##    message:
##        'STAR - Creating: {output} '
##    run:
##        assert os.path.exists(input.star_index)
##        shell('''
##        STAR \
##        --genomeDir {params.star_index} \
##        --genomeLoad LoadAndKeep \
##        --readFilesIn {input.R1} {input.R2} \
##        --readFilesCommand gunzip -c \
##        --outFileNamePrefix {params.out_prefix} \
##        --outSAMtype BAM Unsorted \
##        ''')
#
##rule se_STAR_mapping:
##    input:
##        R1 = '/scratch/single_end_mapping/trimmed_data/{sample}_se_trim.fastq.gz',
##        #star_index = '/scratch/HorseGeneAnnotation/public/refgen/GCF_002863925.1_EquCab3.0/STAR_INDICES/download.done'
##    params:
##        out_prefix = '/scratch/single_end_mapping/star_map_se/RNASEQ/bam/{sample}_se_',
##        star_index = '/scratch/RNAMapping/HorseGeneAnnotation/public/refgen/GCF_002863925.1_EquCab3.0/STAR_INDICES'
##    output:
##        '/scratch/single_end_mapping/star_map_se/RNASEQ/bam/{sample}_se_Aligned.out.bam'
##    message:
##        'STAR - Creating: {output} '
##    run:
##        #assert os.path.exists(params.star_index)
##        shell('''
##        STAR \
##        --genomeDir {params.star_index} \
##        --genomeLoad LoadAndKeep \
##        --readFilesIn {input.R1} \
##        --readFilesCommand gunzip -c \
##        --outFileNamePrefix {params.out_prefix} \
##        --outSAMtype BAM Unsorted \
##        ''')
#
## ----------------------------------------------------------
##       StringTie
## ----------------------------------------------------------
#
# SORT BAMS
rule pe_sort_bam:
    input:
        bam = S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/bam/{GCF}/paired_end/{sample}_Aligned.out.bam')
    output:
        sortedbam = S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/bam/{GCF}/paired_end/{sample}.sorted.bam')
    shell:
        '''
        samtools view -u {input.bam} | samtools sort - -o {output.sortedbam}
        '''
#
##rule se_sort_bam:
##    input:
##        bam = '/scratch/single_end_mapping/star_map_se/RNASEQ/bam/{sample}_se_Aligned.out.bam'
##    output:
##        sorted_bam = '/scratch/single_end_mapping/star_map_se/RNASEQ/bam/{sample}_se.sorted.bam'
##    shell:
##        '''
##        samtools view -u {input.bam} | samtools sort -o {output.sorted_bam}
##        '''
#
# RUN STRINGTIE
rule pe_run_stringtie:
    input:
        bam = S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/bam/{GCF}/paired_end/{sample}.sorted.bam'),
        ref_gff = S3.remote('HorseGeneAnnotation/public/refgen/{GCF}/{GCF}_genomic.nice.gff.gz',keep_local=True)
    output:
        gff = 'HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/{sample}.gff'
    log:
	'HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/{sample}.log'
    run:
        shell('''
        stringtie \
        {input.bam} \
        -G {input.ref_gff} \
        -o {output.gff}
        ''')
#
##rule se_run_stringtie:
##    input:
##        bam = '/scratch/single_end_mapping/star_map_se/RNASEQ/bam/{sample}_se.sorted.bam',
##        gff = expand('/scratch/RNAMapping/HorseGeneAnnotation/public/refgen/{GCF}/{GCF}_genomic.nice.gff.gz', GCF=config['GCF'])
##    output:
##        gff = f"/scratch/single_end_mapping/refgen/{config['GCF']}/GFF/{{sample}}.gff"
##    run:
##        shell('''
##        stringtie \
##        {input.bam} \
##        -G {input.gff} \
##        -o {output.gff}
##        ''')
#
rule create_merged_GFF_file:
    input:
#        gffs = S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/{sample}.gff',sample=SAMPLES,GCF=REF_GFF))
        gffs = expand('HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/{sample}.gff',sample=SAMPLES,GCF=REF_GFF)
    output:
        gff_list = 'HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/all_GFFs_list.txt',
    log:
	'HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/all_GFFs_list.log'
    run:
        with open(output.gff_list,'w') as OUT:
            print('\n'.join(input.gffs),file=OUT)
#
#
# STRINGTIE ON MERGE
rule pe_stringtie_merge:
    input:
#        gffs = S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/{sample}.gff',sample=SAMPLES,GCF=REF_GFF)),
        gffs = expand('HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/{sample}.gff',sample=SAMPLES,GCF=REF_GFF),
	ref_gff = S3.remote('HorseGeneAnnotation/public/refgen/{GCF}/{GCF}_genomic.nice.gff.gz'),
        gff_list = 'HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/all_GFFs_list.txt'
    output:
        merged_gff = 'HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/Merged.gff'
    log:
	'HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/Merged.log'
    run:
        shell('''
            stringtie \
            --merge -p 10 -o {output.merged_gff} \
            -G {input.ref_gff} \
            {input.gff_list}
        ''')
#
##rule se_stringtie_merge:
##    input:
##        gffs = expand(f"/scratch/single_end_mapping/refgen/{config['GCF']}/GFF/{{sample}}.gff", sample=SE_SAMPLES),
##        ref_gff = expand('/scratch/RNAMapping/HorseGeneAnnotation/public/refgen/{GCF}/{GCF}_genomic.nice.gff.gz', GCF=config['NCBI_GCF'])
##    output:
##        merged_gff = f"/scratch/single_end_mapping/refgen/{config['GCF']}/GFF/merged.gff"
##    run:
##        with open('all_GFFs_list.txt','w') as OUT:
##            print('\n'.join(input.gffs),file=OUT)
##        shell('''
##        stringtie \
##        --merge -p 10 -o {output.merged_gff} \
##        -G {input.ref_gff} \
##        all_GFFs_list.txt
##        ''')
#
# STRINGTIE RECALC ON MERGE
rule pe_stringtie_recalculate_on_merged:
    input:
        #bam = expand('HorseGeneAnnotation/private/sequence/RNASEQ/bam/{GCF}/paired_end/{sample}.sorted.bam',sample=SAMPLES,GCF=REF_GFF),
        #merged_gff = expand('HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/Merged.gff',GCF=REF_GFF)
        bam = 'HorseGeneAnnotation/private/sequence/RNASEQ/bam/{GCF}/paired_end/{sample}.sorted.bam',
        merged_gff = 'HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/Merged.gff'
    output:
        countsdir = directory('HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/Merged/{sample}/counts'),
        gff = S3.remote('HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/Merged/{sample}/{sample}.gff',keep_local=True)
    run:
        shell('''
        stringtie \
        -e \
        -b {output.countsdir}\
        -G {input.merged_gff} \
        -o {output.gff} \
        {input.bam}
        ''')
#
##rule se_stringtie_recalc_on_merged:
##    input:
##        bam = '/scratch/single_end_mapping/star_map_se/RNASEQ/bam/{sample}_se.sorted.bam',
##        merged_gff = f"/scratch/single_end_mapping/refgen/{config['GCF']}/GFF/merged.gff"
##    output:
##        counts_dir = directory(f"/scratch/single_end_mapping/refgen/{config['GCF']}/GFF/Merged/{{sample}}/counts"),
##        gff = f"/scratch/single_end_mapping/refgen/{config['GCF']}/GFF/Merged/{{sample}}/{{sample}}.gff"
##    run:
##        shell('''
##        stringtie \
##        -e \
##        -b {output.counts_dir} \
##        -G {input.merged_gff} \
##        -o {output.gff} \
##        {input.bam}
##        ''')
#
#
## ----------------------------------------------------------
##       Make FPKM tables
## ----------------------------------------------------------
#
##rule download_counts:
##    input:
##        'HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/Merged/{sample}/counts/download.done'
##
##
##rule stringtie_counts:
##    input:
##        'HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/Merged/{sample}/counts/t_data.ctab'
##    output:
##        touch('HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/Merged/{sample}/counts/download.done')
#
#
rule pe_make_FPKM_tables:
    input:
        #counts_dl = expand('HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/Merged/{sample}/counts/download.done',GCF=REF_GFF,sample=SAMPLES)
        counts = expand('HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/Merged/{sample}/counts',GCF=REF_GFF,sample=SAMPLES)
    output:
        transcript_fpkm = 'HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/Merged/transcript_fpkm.tsv',
        gene_fpkm = 'HorseGeneAnnotation/public/refgen/{GCF}/paired_end/GFF/Merged/gene_fpkm.tsv'
    run:
        import pandas as pd
        import numpy  as np
        dfs = []
        ctabs = [os.path.join(c,'t_data.ctab') for c in input.counts]
        for sample,f in zip(SAMPLES,ctabs):
            df = pd.read_table(f)
            df['sample'] = sample
            dfs.append(df)
        df = pd.concat(dfs)
        by_transcript = pd.pivot_table(df,index='t_name',columns='sample',values='FPKM')
        by_transcript.to_csv(output.transcript_fpkm,sep='\t')
        by_gene = pd.pivot_table(df,index='gene_name',columns='sample',values='FPKM',aggfunc=np.mean)
        by_gene.to_csv(output.gene_fpkm,sep='\t')

#rule se_make_FPKM_tables:
#    input:
#        counts = expand('/scratch/single_end_mapping/refgen/{GCF}/GFF/Merged/{sample}/counts/t_data.ctab', GCF=config['GCF'], sample=SE_SAMPLES)
#    output:
#        transcript_fpkm = f"/scratch/single_end_mapping/refgen/{config['GCF']}/GFF/Merged/transcript_fpkm.tsv",
#        gene_fpkm = f"/scratch/single_end_mapping/refgen/{config['GCF']}/GFF/Merged/gene_fpkm.tsv"
#    run:
#        import pandas as pd
#        import numpy as np
#        dfs = []
#        for sample,f in zip(SE_SAMPLES,input):
#            df = pd.read_table(f)
#            df['sample'] = sample
#            dfs.append(df)
#        df = pd.concat(dfs)
#        by_transcript = pd.pivot_table(df,index='t_name',columns='sample',values='FPKM')
#        by_transcript.to_csv(output.transcript_fpkm,sep='\t')
#        by_gene = pd.pivot_table(df,index='gene_name',columns='sample',values='FPKM',aggfunc=np.mean)
#        by_gene.to_csv(output.gene_fpkm,sep='\t')

# ----------------------------------------------------------
#       Make TPM matrix
# ----------------------------------------------------------

# need index

rule SALMON_mapping:
    input:
        R1 = 'HorseGeneAnnotation/private/sequence/trimmed/{sample}_trim1.fastq.gz',
        R2 = 'HorseGeneAnnotation/private/sequence/trimmed/{sample}_trim2.fastq.gz',
        refgen = '/scratch/NovaseqComp/novaseq/HorseGeneAnnotation/public/refgen/GCF_002863925.1_EquCab3.0/SALMON_INDEX'
    output:
#        directory('HorseGeneAnnotation/public/refgen/{GCF}/paired_end/SalmonQuant/{sample}'),
        'HorseGeneAnnotation/public/refgen/{GCF}/paired_end/SalmonQuant/{sample}/quant.sf'
    params:
        out_dir = 'HorseGeneAnnotation/public/refgen/{GCF}/paired_end/SalmonQuant/{sample}'
    run:
#        out_dir = os.path.dirname(output)
        shell('''
            salmon quant \
            -i {input.refgen} \
            -l A \
            -1 {input.R1} \
            -2 {input.R2} \
            --validateMappings \
            -o {params.out_dir}
        ''')

rule make_TPM_matrix:
    input:
        quant_files = expand('HorseGeneAnnotation/public/refgen/{GCF}/paired_end/SalmonQuant/{sample}/quant.sf',sample=SAMPLES,GCF=REF_GFF)
    output:
        tpm_matrix = 'HorseGeneAnnotation/public/refgen/{GCF}/paired_end/SalmonQuant/novaseq_TPM.tsv',
#        fat_matrix = 'HorseGeneAnnotation/public/refgen/{GCF}/paired_end/SalmonQuant/Fat_TPM.tsv',
#        muscle_matrix = 'HorseGeneAnnotation/public/refgen/{GCF}/paired_end/SalmonQuant/Muscle_TPM.tsv'
    run:
        tpm = None
        for n in input.quant_files:                   
            name = n.split('/')[-2]
            print(name)
            #x = pd.read_table(f"data/SalmonQuant/{n}/quant.sf")                                         
            x = pd.read_table(n)
            x = x.loc[:,['Name','TPM']]            
            x.rename(columns={'Name':'gene','TPM':name},inplace=True)                                      
            x.set_index('gene',inplace=True)       
            if tpm is None:                        
                tpm = x                            
            else:                                  
                tpm = tpm.join(x,lsuffix='gene') 
        # Create the whole matrix
        tpm.to_csv(output.tpm_matrix,sep='\t')
        # Create the Fat & Muscle Matrix
#        tpm.loc[:,[x.endswith('F') for x in tpm.columns]].to_csv(output.fat_matrix,sep='\t')
#        tpm.loc[:,[x.endswith('M') for x in tpm.columns]].to_csv(output.muscle_matrix,sep='\t')
