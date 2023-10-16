import subprocess
import os
import shutil

if os.path.isfile("gtdbtk/classify/gtdbtk.bac120.summary.tsv") and os.path.isfile("gtdbtk/classify/gtdbtk.ar53.summary.tsv"):
    subprocess.Popen(
        """
    	DRAM.py annotate \
    	-i spades_assembly/scaffolds.fasta \
    	--threads %s \
    	--gtdb_taxonomy gtdbtk/classify/gtdbtk.bac120.summary.tsv \
    	--gtdb_taxonomy gtdbtk/classify/gtdbtk.ar53.summary.tsv \
    	-o dram/annotation
        """ % (snakemake.threads), shell=True).wait()

elif os.path.isfile("gtdbtk/classify/gtdbtk.bac120.summary.tsv") and not os.path.isfile("gtdbtk/classify/gtdbtk.ar53.summary.tsv"):
    subprocess.Popen(
        """
    	DRAM.py annotate \
    	-i spades_assembly/scaffolds.fasta \
    	--threads %s \
    	--gtdb_taxonomy gtdbtk/classify/gtdbtk.bac120.summary.tsv \
    	-o dram/annotation
        """ % (snakemake.threads), shell=True).wait()

elif os.path.isfile("gtdbtk/classify/gtdbtk.ar53.summary.tsv") and not os.path.isfile("gtdbtk/classify/gtdbtk.bac120.summary.tsv"):
    subprocess.Popen(
        """
    	DRAM.py annotate \
    	-i spades_assembly/scaffolds.fasta \
    	--threads %s \
    	--gtdb_taxonomy gtdbtk/classify/gtdbtk.ar53.summary.tsv \
    	-o dram/annotation
        """ % (snakemake.threads), shell=True).wait()




if os.path.isfile("dram/annotation/trnas.tsv") and os.path.isfile("dram/annotation/rrnas.tsv"):
    subprocess.Popen(
        """
        DRAM.py distill \
        -i dram/annotation/annotations.tsv \
        -o dram/genome_summaries \
        --trna_path dram/annotation/trnas.tsv \
        --rrna_path dram/annotation/rrnas.tsv && 
        touch dram/done
        """, shell=True).wait()

elif os.path.isfile("dram/annotation/trnas.tsv") and not os.path.isfile("dram/annotation/rrnas.tsv"):
    subprocess.Popen(
        """
        DRAM.py distill \
        -i dram/annotation/annotations.tsv \
        -o dram/genome_summaries \
        --trna_path dram/annotation/trnas.tsv && 
        touch dram/done
        """, shell=True).wait()

elif os.path.isfile("dram/annotation/rrnas.tsv") and not os.path.isfile("dram/annotation/trnas.tsv"):
    subprocess.Popen(
        """
        DRAM.py distill \
        -i dram/annotation/annotations.tsv \
        -o dram/genome_summaries \
        --rrna_path dram/annotation/rrnas.tsv && 
        touch dram/done
        """, shell=True).wait()





