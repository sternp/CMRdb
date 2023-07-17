import subprocess
import os
import shutil

short_reads_1 = snakemake.config["short_reads_1"]
short_reads_2 = snakemake.config["short_reads_2"]
long_reads = snakemake.config["long_reads"]
output_dir_path = os.path.join(snakemake.config["output"], "qc")
max_memory = str(snakemake.config["max_memory"]) + 'g' # kneaddata takes bytes as input. change to GBs
sequencer_source = snakemake.config["sequencer_source"]
skip_qc = snakemake.config["skip_qc"]
trimmomatic_exec = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../bin/Trimmomatic-0.39/")


print('Short R1 = ' + short_reads_1)
print('Short R2 = ' + short_reads_2)
print('Long = ' + long_reads)
shutil.rmtree("qc/")


if skip_qc == True and long_reads == "none":
    subprocess.Popen(
        """
        echo '--skip-qc specified, skipping read QC steps...' &&
        mkdir qc
        mkdir qc/clean_reads &&
        mkdir qc/clean_reads/R1 &&
        mkdir qc/clean_reads/R2 &&
        cp %s qc/clean_reads/R1 &&
        cp %s qc/clean_reads/R2 &&
        touch qc/done
        """ %
        (short_reads_1, short_reads_2), shell=True).wait()

elif skip_qc == True and long_reads != "none":
    subprocess.Popen(
        """
        echo '--skip-qc specified, skipping read QC steps...' &&
        mkdir qc
        mkdir qc/clean_reads &&
        mkdir qc/clean_reads/R1 &&
        mkdir qc/clean_reads/R2 &&
        mkdir qc/clean_reads/long &&
        cp %s qc/clean_reads/R1 &&
        cp %s qc/clean_reads/R2 &&
        cp %s qc/clean_reads/long &&
        touch qc/done
        """
        , shell=True).wait()

elif short_reads_1 != "none" and short_reads_2 != "none" and long_reads == "none" and skip_qc != True:
    subprocess.Popen(
        """
        kneaddata \
        -i1 %s \
        -i2 %s \
        -o %s \
        -t %s \
        -p %s \
        --max-memory %s \
        --sequencer-source %s \
        --trimmomatic %s \
        --run-fastqc-end \
        --reorder \
        --remove-intermediate-output \
        --bypass-trf \
		--run-trim-repetitive \
        --decontaminate-pairs strict &&

        mkdir qc/clean_reads &&
        mkdir qc/clean_reads/R1/ &&
        mkdir qc/clean_reads/R2/ &&

        mv qc/*trimmed.1.fastq qc/clean_reads/R1 &&
        mv qc/*trimmed.2.fastq qc/clean_reads/R2 &&

        pigz qc/clean_reads/R1/* &&
        pigz qc/clean_reads/R2/* &&

        touch qc/done
        """ %
        (short_reads_1, short_reads_2, output_dir_path, snakemake.threads, snakemake.threads, max_memory, sequencer_source, trimmomatic_exec), shell=True).wait()

elif short_reads_1 != "none" and short_reads_2 != "none" and long_reads != "none" and skip_qc != True:
    subprocess.Popen(
        """
        kneaddata \
        -i1 %s \
        -i2 %s \
        -o %s \
        -t %s \
        -p %s \
        --max-memory %s \
        --sequencer-source %s \
        --trimmomatic %s \
        --run-fastqc-end \
        --reorder \
        --remove-intermediate-output \
        --bypass-trf \
		--run-trim-repetitive \
        --decontaminate-pairs strict &&

        mkdir qc/clean_reads &&
        mkdir qc/clean_reads/R1/ &&
        mkdir qc/clean_reads/R2/ &&

        mv qc/*trimmed.1.fastq qc/clean_reads/R1 &&
        mv qc/*trimmed.2.fastq qc/clean_reads/R2 &&

        pigz qc/clean_reads/R1/* &&
        pigz qc/clean_reads/R2/* &&

        """ %
        (short_reads_1, short_reads_2, output_dir_path, snakemake.threads, snakemake.threads, max_memory, sequencer_source, trimmomatic_exec), shell=True).wait()

    subprocess.Popen(
        """
        mkdir qc/clean_reads/long/ &&
        cp %s qc/clean_reads/long/ &&

        touch qc/done
        """ %
        (long_reads), shell=True).wait()
