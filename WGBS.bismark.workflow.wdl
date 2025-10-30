version 1.0

workflow_WGBS
{
    input {
        String sample
        File raw_fq1
        File raw_fq2
        File target_bed
        String output_dir

        # options
        # Mapping efficiency
        Int map_eff = 30
        # depth
        Int depth = 5
        File bowtie2_dir
        String bismark_param = "-bowtie2 -L 32 -N 1 -non_directional"
        File target_bin_dir
    }

    String docker_img = 'docker_image'

    # Quality control of fastq
    call fastq_qc {
        input:
            sample       = sample,
            fq1          = raw_fq1,
            fq2          = raw_fq2,
            out_dir      = "${output_dir}/${sample}/Fastqc",
            cpu          = 1,
            mem          = 1,
            docker_image = docker_img
    }

    # Trimmed the adapters and low-quality bases by trim_galore with default parameters.
    call preprocess {
        input:
            sample       = sample,
            fq1          = raw_fq1,
            fq2          = raw_fq2,
            out_dir      = "${output_dir}/${sample}/Trim_galore",
            cpu          = 1,
            mem          = 1,
            docker_image = docker_img
    }

    # Created Bisulfite Genome, aligned the reads to hg19/hg38 using bismark
    # Deduplicate && report 
    call bismark_cpg {
        input:
            sample              = sample,
            filter_fq1          = preprocess.filter_fq1,
            filter_fq2          = preprocess.filter_fq2,
            bowtie2_dir         = bowtie2_dir,
            bismark_param       = ${bismark_param},
            out_dir             = "${output_dir}/${sample}/Bismark",
            cpu                 = 1,
            mem                 = 1,
            docker_image        = docker_img
    }

    # Extract the methylation information of target region && statistic
    call target_cpg {
        input:
            sample         = sample,
            coverage       = bismark_cpg.cov,
            bam            = bismark_cpg.bam,
            report         = bismark_cpg.report,
            target_bin_dir = target_bin_dir,
            bed            = ${target_bed},
            dp             = ${depth},
            out_dir        = "${output_dir}/${sample}/Bismark",
            cpu            = 1,
            mem            = 1,
            docker_image   = docker_img
    }

    output {
        # fastqc
        File raw_fq1_qc = fastq_qc.fq1_qc
        File raw_fq2_qc = fastq_qc.fq2_qc

        # preprocess of fastq
        File clean_fq1 = preprocess.filter_fq1
        File clean_fq2 = preprocess.filter_fq2
        File trim_fq1_report = preprocess.fq1_trim_report
        File trim_fq2_report = preprocess.fq2_trim_report

        # bismark genome, alignment, deduplicate and report
        File all_cov_file = bismark_cpg.cov
        File all_bam_file = bismark_cpg.bam
        File all_report_file = bismark_cpg.report

        # target
        File target_cov_file = target_cpg.target_cov
        File target_stats_file = target_cpg.target_stats
    }
}


task fastq_qc {
    # qc of fastq
    input {
        String sample
        File fq1
        File fq2
        String out_dir

        # Runtime environment setting
        Int cpu  = 1
        Int mem  = 1
        String docker_image = docker_img
    }

    command <<<
        mkdir -p ~{out_dir}
        fastqc ~{fq1} ~{fq2} --outdir ~{out_dir}
        out_prefix1=$(basename ~{fq1} | sed -r 's/.fastq.gz|.fq.gz//g')
        out_prefix2=$(basename ~{fq2} | sed -r 's/.fastq.gz|.fq.gz//g')
    >>>

    output {
        File fq1_qc = "~{out_dir}/${out_prefix1}_fastqc.html"
        File fq2_qc = "~{out_dir}/${out_prefix2}_fastqc.html"
    }

    runtime {
        req_cpu: cpu
        req_memory: "${mem}G"
        docker_url: docker_image
    }
}


task preprocess {
    # Trimmed the adapters and low-quality bases
    input {
        String sample
        File fq1
        File fq2
        String out_dir
        
        # Runtime environment setting
        Int cpu = 1
        Int mem = 1
        String docker_image = docker_img
    }

    command {
        mkdir -p ${out_dir}
        trim_galore --illumina --gzip \
            --basename ${sample} \
            --output_dir ${out_dir} \
            --cores ${cpu} \
            --paired ${fq1} ${fq2} 2>${out_dir}/${sample}.trim_galore.log
    }

    output {
        File filter_fq1 = "${out_dir}/${sample}_val_1.fq.gz"
        File filter_fq2 = "${out_dir}/${sample}_val_2.fq.gz"
        File fq1_trim_report = "${out_dir}/${fq1}_trimming_report.txt"
        File fq2_trim_report = "${out_dir}/${fq2}_trimming_report.txt"
    }

    runtime {
        req_cpu: cpu
        req_memory: "${mem}G"
        docker_url: docker_image
    }
}


task bismark_cpg {
    # Created Bisulfite Genome, aligned the reads to hg19/hg38 using bismark
    # Deduplicate && report 
    input {
        String sample
        File filter_fq1
        File filter_fq2
        Int map_eff = 30 
        String bismark_param
        File bowtie2_dir
        String out_dir

        # Runtime environment setting
        Int cpu = 1
        Int mem = 1
        String docker_image  
    }

    comand <<<
        mkdir -p ~{out_dir}/genome

        # After bismark_genome_preparation construct the index, it could be used for further analysis. 
        # ~{bowtie2_dir}: define the bowtie2 install directory with yourself

        bismark_genome_preparation --path_to_aligner ~{bowtie2_dir} --verbose ~{out_dir}/genome 2>~{out_dir}/genome/~{sample}.bismark_genome_preparation.log

        bismark ~{bismark_param} --genome_folder ~{out_dir}/genome -1 ~{filter_fq1} -2 ~{filter_fq2} --output_dir ~{out_dir} --temp_dir ~{out_dir}
        
        grep "Mapping efficiency" ~{out_dir}/~{sample}_val_1_bismark_bt2_PE_report.txt| sed 's/\%//g'|awk '{if($NF < '${map_eff}') {print "Mapping efficiency < '${map_eff}'%, it suggest to adjust the '"-L & -N"' option of bowtie2 and aligned again";exit 0;} else {}}'

        # Deduplicate && report 
        deduplicate_bismark --paired --bam ~{out_dir}/~{sample}_val_1_bismark_bt2_pe.bam --output_dir ~{out_dir}

        bismark_methylation_extractor --gzip --no_overlap --bedGraph --comprehensive --cytosine_report --parallel ~{cpu} --genome_folder ~{out_dir}/genome --output_dir ~{out_dir} ~{out_dir}/~{sample}_val_1_bismark_bt2_pe.deduplicated.bam 2>~{out_dir}/~{sample}.bismark_methylation_extractor.log

        cd ~{out_dir} && bismark2report --dir ~{out_dir}
        bismark2summary --basename ~{sample}.bismark_summary_report
    >>>

    output {
        File cov = "~{out_dir}/~{sample}_val_1_bismark_bt2_pe.deduplicated.bismark.cov"
        File bam = "~{out_dir}/~{sample}_val_1_bismark_bt2_pe.deduplicated.bam"
        File report = "~{out_dir}/~{sample}_val_1_bismark_bt2_PE_report.txt"
    }

    runtime {
        req_cpu: cpu
        req_memory: "${mem}G"
        docker_url: docker_image
    }
}


task target_cpg {
    # Extract the methylation information of target region && statistic
    input {
        String sample
        File coverage
        File bam
        File report
        File bed
        File target_bin_dir
        Int dp = 5
        String out_dir

        # Runtime environment setting
        Int cpu = 1 
        Int mem = 1
        String docker_image = docker_img    
    }

    command <<<
        mkdir -p ~{out_dir}/Statistic

        # target methylation cov
        python ~{target_bin_dir}/extract.target.met.cov.py --sample c --dp ~{dp} --cov ~{coverage} --bed ~{bed} --outDir ~{out_dir}

        # statistic
        mapRatio=$(grep "Mapping efficiency" ~{report} | awk '{print $NF}')

        samtools view -bh -L ~{bed} ~{bam} > ~{out_dir}/~{sample}_val_1_bismark_bt2_pe.deduplicated.target.bam

        bedtools coverage -a ~{bed} -b ~{bam} > ~{out_dir}/Statistic/~{sample}.deduplicated.bam.target.bed.cov.txt

        cov=$(cat ~{out_dir}/Statistic/~{sample}.deduplicated.bam.target.bed.cov.txt | awk '{sum1+=$5; sum2+=$6} END {print sum1*100/sum2}')
        meanDep=$(cat ~{out_dir}/Statistic/~{sample}.deduplicated.bam.target.bed.cov.txt | awk '{sum1+=$4; sum2+=$6} END {print sum1*100/sum2}')

        samtools sort ~{out_dir}/~{sample}_val_1_bismark_bt2_pe.deduplicated.target.bam -o ~{out_dir}/~{sample}_val_1_bismark_bt2_pe.deduplicated.target.sort.bam
        samtools index ~{out_dir}/~{sample}_val_1_bismark_bt2_pe.deduplicated.target.sort.bam
        samtools depth ~{out_dir}/~{sample}_val_1_bismark_bt2_pe.deduplicated.target.sort.bam > ~{out_dir}/Statistic/~{sample}.deduplicated.bam.target.bed.depth.txt

        num=$(wc -l ~{out_dir}/Statistic/~{sample}.deduplicated.bam.target.bed.depth.txt|awk '{print $(NF-1)}')

        # depth 1x,5x,10x,30x
        cov1Dep=$(awk -F "\t" 'BEGAIN{b=0}{if($3>=1){b+=$3}}END{print b/"'$num'"}' ~{out_dir}/Statistic/~{sample}.deduplicated.bam.target.bed.depth.txt)
        cov5Dep=$(awk -F "\t" 'BEGAIN{b=0}{if($3>=5){b+=$3}}END{print b/"'$num'"}' ~{out_dir}/Statistic/~{sample}.deduplicated.bam.target.bed.depth.txt)
        cov10Dep=$(awk -F "\t" 'BEGAIN{b=0}{if($3>=10){b+=$3}}END{print b/"'$num'"}' ~{out_dir}/Statistic/~{sample}.deduplicated.bam.target.bed.depth.txt)
        cov30Dep=$(awk -F "\t" 'BEGAIN{b=0}{if($3>=30){b+=$3}}END{print b/"'$num'"}' ~{out_dir}/Statistic/~{sample}.deduplicated.bam.target.bed.depth.txt)

        echo -e "Sample\tMapping Ratio\tOn-target coverage (%)\tOn-target average depth (x)\tOn-target 1x\tOn-target 5x\tOn-target 10x\tOn-target 30x" > ~{out_dir}/Statistic/~{sample}_val_1_bismark_bt2_pe.deduplicated.target.statistic.txt
        echo -e ~{sample}"\t$mapRatio\t$cov\t$meanDep\t$cov1Dep\t$cov5Dep\t$cov10Dep\t$cov30Dep" >> ~{out_dir}/Statistic/~{sample}_val_1_bismark_bt2_pe.deduplicated.target.statistic.txt
    >>>

    output {
        File target_cov = "~{out_dir}/~{sample}.target.bismark.filter.cov.gz"
        File target_stats = "~{out_dir}/Statistic/~{sample}_val_1_bismark_bt2_pe.deduplicated.target.statistic.txt"
    }

    runtime {
        req_cpu: cpu
        req_memory: "${mem}G"
        docker_url: docker_image
    }
}