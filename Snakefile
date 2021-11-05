#!/usr/bin/env python3

import pathlib2
import pandas
import os
import peppy

#############
# FUNCTIONS #
#############

def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))

def find_read_files(read_dir):
#Make list of files
    path_generator = os.walk(read_dir, followlinks = True)
    my_files = list((dirpath, filenames)
        for (dirpath, dirname, filenames)
        in path_generator)
#Make new dictionary & populate with files (flowcell = key)
    my_fastq_files = {}
    for dirpath, filenames in my_files:
        for filename in filenames:
            if filename.endswith('.fastq.gz'):
                my_flowcell = pathlib2.Path(dirpath).name
                my_fastq = str(pathlib2.Path(dirpath,filename))
                if my_flowcell in my_fastq_files:
                    my_fastq_files[my_flowcell].append(my_fastq)
                else:
                    my_fastq_files[my_flowcell]= []
                    my_fastq_files[my_flowcell].append(my_fastq)
    return(my_fastq_files)

def sample_name_to_fastq(wildcards):
    sample_row = sample_key[sample_key['Sample_name'] == wildcards.sample]
    sample_id = sample_row.iloc[-1]['OGF_sample_ID']
    sample_flowcell = sample_row.iloc[-1]['Flow_cell']
    sample_all_fastq = [x for x in all_fastq[sample_flowcell]
                        if '-{}-'.format(sample_id) in x]
    sample_r1 = sorted(list(x for x in sample_all_fastq
                            if '_R1_' in os.path.basename(x)))
    sample_r2 = sorted(list(x for x in sample_all_fastq
                            if '_R2_' in os.path.basename(x)))
    return({'r1': sample_r1, 'r2': sample_r2})

###########
# GLOBALS #
###########

##this parses the config & sample key files into an object named pep
pepfile: 'data/config.yaml'
##can now use this to generate list of all samples
all_dna_samples=pep.sample_table["sample_name"]


read_dir = 'data/rna_reads'
sample_key_file = 'data/sample_key.csv'

bbduk_adapters = '/adapters.fa'
bbduk_ref = '/phix174_ill.ref.fa.gz'

#containers
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'

#########
# SETUP #
#########

##https://benlangmead.github.io/aws-indexes/k2 for kraken db viral (5/17/2021)

# generate name to filename dictionary
all_fastq = find_read_files(read_dir)

sample_key = pandas.read_csv(sample_key_file)
all_rna_samples = sorted(set(sample_key['Sample_name']))

#########
# RULES #
#########

rule target:
    input:
        expand('output/kraken_rna/{sample}/{sample}_kraken_report.txt', sample=all_rna_samples),
        expand('output/kraken_dna/{sample}/{sample}_kraken_report.txt', sample=all_dna_samples),
        'output/kraken_dna/genome/genome_kraken_report.txt'

##########################
## Genome DNAseq Kraken ##
##########################

rule kraken_genome:
    input:
        r1 = 'output/bbduk_trim_dna/genome/genome_trimr1.fq.gz',
        r2 = 'output/bbduk_trim_dna/genome/genome_trimr2.fq.gz',
        db = 'kraken_db/viral'
    output:
        out = 'output/kraken_dna/genome/genome_kraken_out.txt',
        report = 'output/kraken_dna/genome/genome_kraken_report.txt'
    log:
        'output/logs/kraken_out_genome.log'
    threads:
        20
    shell:
        'bin/kraken_v2.1.2/kraken2 '
        '--threads {threads} '
        '--db {input.db} '
        '--paired '
        '--output {output.out} '
        '--report {output.report} '
        '--use-names '
        '{input.r1} {input.r2} '
        '&> {log}'

##trim and decontaminate DNA reads to map onto genome
rule bbduk_trim_genome:
    input:
        filr1 = 'output/bbduk_trim_dna/genome/genome_filr1.fq.gz',
        filr2 = 'output/bbduk_trim_dna/genome/genome_filr2.fq.gz'
    output:
        trimr1 = 'output/bbduk_trim_dna/genome/genome_trimr1.fq.gz',
        trimr2 = 'output/bbduk_trim_dna/genome/genome_trimr2.fq.gz',
        t_stats = 'output/bbduk_trim_dna/genome/genome/trim-stats.txt'
    log:
        trim = 'output/logs/bbduk_trim_dna/genome_trim.log'
    params:
        adapters = bbduk_adapters
    threads:
        20
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in1={input.filr1} '
        'in2={input.filr2} '
        'int=f '
        'out1={output.trimr1} '
        'out2={output.trimr2} '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.t_stats} '
        '2> {log.trim} '

rule bbduk_filter_genome:  
    input:
        r1 = 'data/dna_reads/genome_r1.fastq.gz',
        r2 = 'data/dna_reads/genome_r2.fastq.gz'
    output:
        filr1 = 'output/bbduk_trim_dna/genome/genome_filr1.fq.gz',
        filr2 = 'output/bbduk_trim_dna/genome/genome_filr2.fq.gz',
        f_stats = 'output/bbduk_trim_dna/genome/genome_filter-stats.txt'
    log:
        filter = 'output/logs/bbduk_trim_dna/genome_filter.log'
    params:
        ref = bbduk_ref
    threads:
        20
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in1={input.r1} '
        'in2={input.r2} '
        'out1={output.filr1} '
        'out2={output.filr2} '
        'ref={params.ref} '
        'hdist=1 '
        'stats={output.f_stats} '       
        '2> {log.filter} ' 

#######################
## Pop DNAseq Kraken ##
#######################

rule kraken_pop:
    input:
        reads = 'output/bbduk_trim_dna/{sample}/{sample}_trim.fq.gz',
        db = 'kraken_db/viral'
    output:
        out = 'output/kraken_dna/{sample}/{sample}_kraken_out.txt',
        report = 'output/kraken_dna/{sample}/{sample}_kraken_report.txt'
    log:
        'output/logs/kraken_out_{sample}.log'
    threads:
        20
    shell:
        'bin/kraken_v2.1.2/kraken2 '
        '--threads {threads} '
        '--db {input.db} '
        '--output {output.out} '
        '--report {output.report} '
        '--use-names '
        '{input.reads} '
        '&> {log}'


##trim and decontaminate DNA reads to map onto genome
rule bbduk_trim_dna:
    input:
        fastq = 'output/bbduk_trim_dna/{sample}/{sample}_filr.fq.gz'
    output:
        trim_fastq = 'output/bbduk_trim_dna/{sample}/{sample}_trim.fq.gz',
        t_stats = 'output/bbduk_trim_dna/{sample}/trim-stats.txt'
    log:
        trim = 'output/logs/bbduk_trim_dna/{sample}_trim.log'
    params:
        adapters = bbduk_adapters
    threads:
        20
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in={input.fastq} '
        'int=t '
        'out={output.trim_fastq} '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.t_stats} '
        '2> {log.trim} '

rule bbduk_filter_dna:  
    input:
        fastq = 'data/dna_reads/{sample}.fastq'
    output:
        fastq = 'output/bbduk_trim_dna/{sample}/{sample}_filr.fq.gz',
        f_stats = 'output/bbduk_trim_dna/{sample}/filter-stats.txt'
    log:
        filter = 'output/logs/bbduk_trim_dna/{sample}_filter.log'
    params:
        ref = bbduk_ref
    threads:
        20
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in={input.fastq} '
        'int=t '
        'out={output.fastq} '
        'ref={params.ref} '
        'hdist=1 '
        'stats={output.f_stats} '       
        '2> {log.filter} ' 

###################
## RNAseq Kraken ##
###################

rule kraken:
    input:
        r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
        r2 = 'output/bbduk_trim/{sample}_r2.fq.gz',
        db = 'kraken_db/viral'
    output:
        out = 'output/kraken_rna/{sample}/{sample}_kraken_out.txt',
        report = 'output/kraken_rna/{sample}/{sample}_kraken_report.txt'
    log:
        'output/logs/kraken_out_{sample}.log'
    threads:
        20
    shell:
        'bin/kraken_v2.1.2/kraken2 '
        '--threads {threads} '
        '--db {input.db} '
        '--paired '
        '--output {output.out} '
        '--report {output.report} '
        '--use-names '
        '{input.r1} {input.r2} '
        '&> {log}'

rule bbduk_trim:
    input:
        r1 = 'output/joined/{sample}_r1.fq.gz',
        r2 = 'output/joined/{sample}_r2.fq.gz'
    output:
        r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
        r2 = 'output/bbduk_trim/{sample}_r2.fq.gz'
    params:
        adapters = bbduk_adapters
    log:
        'output/logs/bbduk_trim/{sample}.log'
    threads:
        20
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.r1} '
        'out2={output.r2} '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=15 '
        '&> {log}'

rule cat_reads:
    input:
        unpack(sample_name_to_fastq)
    output: 
        r1 = temp('output/joined/{sample}_r1.fq.gz'),
        r2 = temp('output/joined/{sample}_r2.fq.gz')
    threads:
        1
    shell:
        'cat {input.r1} > {output.r1} & '
        'cat {input.r2} > {output.r2} & '
        'wait'