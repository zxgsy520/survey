
import os.path
from collections import OrderedDict

ROOT = "/nextomics/Pipeline/Survey/v1.1.1/"
SCRIPTS = os.path.join(ROOT, "scripts")
TEMPLATES = os.path.join(ROOT, "templates") 
BIN = os.path.join(ROOT, "survey")
TAXONOMY = os.path.join(TEMPLATES, "species.taxonomy")

#PYNTHON_bin
PYTHON_BIN="/export/software/Bases/Python/v2.7.15/bin"
#R
RSCRIPT="/export2/master2/xuwj/software/conda/envs/py2.7/bin/Rscript"
#FASTP
FASTP="/export2/master2/zhangxg/script/fastp"
#FASTQC
FASTQC="/nextomics/Software/QC/fastqc/v0.11.8/fastqc"
#KmerFreq
KMERFREQ="/export2/software/Assembly/NGS/SOAPdenovo2-src-r240/KmerFreq_AR"
#Seqkit
SEQKIT="/export2/master2/zhangxg/script/seqkit"
#JELLYFISH
JELLYFISH_BIN="/export2/master2/zhangxg/software/jellyfish/2.3.0"
#SOAPdenovo
SOAPDENOVO_BIN="/export2/software/Assembly/NGS/SOAPdenovo2-src-r240/"
#BWA
BWA_BIN="/export2/software/Bases/bwa/v0.7.15/"
#MINIMAP
MINIMAP_BIN = "/software/minimap2"
#SAMBLASTER
SAMBLASTER_BIN = "/nextomics/Software/Base/samblaster/v0.1.25"
#SAMTOOLS
SAMTOOLS_BIN="/export2/software/Bases/samtools/v1.4/bin/"
#BLAST
BLAST_BIN="/export2/master2/zhangxg/software/blast/2.9.0/bin/"
#NT_DATABASE
NT_DATABASE="/nextomics/Databases/20200330_NT_Database/nt"
#GHOSTSCRIPT
GHOSTSCRIPT="/export2/master2/zhangxg/software/ghostscript/gs"


SEQUENCER = OrderedDict([
    ("RSII", {"minimap2": "-ax map-pb"}
    ),
    ("Sequel", {"minimap2": "-ax map-pb"}
    ),
    ("GridION", {"minimap2": "-ax map-ont"}
    ),
    ("PromethION", {"minimap2": "-ax map-ont"}
    ),
    ("illumina", {"minimap2": "-ax sr"}
    ),
    ("mgi", {"minimap2": "-ax sr"}
    ),
])


# SOFTWARE VERSION

SOFTWARE_VERSION = {
    "fastp":{
        "GETVER": "%s --version 2>&1 |grep 'fastp'" % FASTP,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "0.20.0",
    },
    "fastqc":{
        "GETVER": "%s --version 2>&1 |grep 'FastQC'" % FASTQC,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "0.11.8",
    },
    "blastn":{
        "GETVER": "%s/blastn -version 2>&1| grep 'blastn'" % BLAST_BIN,
        "REGEXP": "\d+\.\d+\.\d+\+",
        "MINVER": "2.7.1+"
    },
    "kmerfreq":{
        "GETVER": "%s --version 2>&1| grep 'Version'" % KMERFREQ,
        "REGEXP": "\d+\.\d+",
        "MINVER": "2.0"
    },
    "jellyfish":{
        "GETVER": "%s/jellyfish --version 2>&1| grep 'jellyfish'" % JELLYFISH_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "2.3.0"
    },
    "soapdenovo":{
        "GETVER": "%s/SOAPdenovo-127mer --version 2>&1| grep 'Version'" % SOAPDENOVO_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "2.04"
    },
    "minimap2": {
        "GETVER": "%s/minimap2 --version" % MINIMAP_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "2.17",
    },
    "samblaster": {
        "GETVER": "%s/samblaster --version 2>&1 |grep 'Version'" % SAMBLASTER_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "0.1.25",
    },
}


# DATABASE CONFIGURE

NT = "/nextomics/Databases/20200330_NT_Database/nt"
NT_TAXON = {
    "plant": NT,
    "animal": NT,
    "fungi": "/nextomics/Databases/nt/20200330/microbe"
}

MIT = "/nextomics/Databases/Refseq/mitochondrion/mitochondrion.genomic.fasta"
PLA = "/nextomics/Databases/Refseq/plastid/plastid.genomic.fasta"
MC_TAXON = {
    "plant": "%s %s" % (MIT, PLA),
    "animal": MIT,
    "fungi": MIT
} 

