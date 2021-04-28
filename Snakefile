# download swiss protein db from: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
DB_FILE = "data/uniprot_sprot.fasta.gz"
CLUSTER = ""
DIAMOND_PATH = "./../diamond-2.0.8"
RAXML_PATH = "raxml-ng"
MUSCLE_PATH = "muscle"


def clean_filename(filename):
    return filename.split("/")[-1].split(".")[0]

DB_FILENAME = clean_filename(DB_FILE)

rule diamond_bin_db:
    input:
        f"{DB_FILE}"
    output:
        f"data/{DB_FILENAME}.dmnd"
    shell:
        "{DIAMOND_PATH}/bin/diamond makedb --in {input} -d {output}"

rule multistep_cluster:
    input:
        f"data/{DB_FILENAME}.dmnd"
    output:
        f"data/{DB_FILENAME}.multistep"
    shell:
        "{DIAMOND_PATH}/bin/diamond cluster --cluster-algo multi-step -d {input} -o {output}"

rule mcl_cluster:
    input:
        f"data/{DB_FILENAME}.dmnd"
    output:
        f"data/{DB_FILENAME}.mcl"
    shell:
        "{DIAMOND_PATH}/bin/diamond cluster --cluster-algo mcl -d {input} -o {output}"

rule separate_cluster:
    input:
        f"{DB_FILE}",
        f"data/{DB_FILENAME}.multistep"
    output:
        f"out/{CLUSTER}.fasta",
        f"out/{CLUSTER}.clean.fasta"
    params:
        cluster_name=CLUSTER
    shell:
        "python3 scripts/clusteread.py w {input} \"{params.cluster_name}\""

rule muscle:
    input:
        f"out/{CLUSTER}.clean.fasta"
    output:
        f"out/{CLUSTER}.msa.fasta"
    shell:
        "{MUSCLE_PATH}/muscle3.8.31_i86darwin64 -in \"{input}\" -out \"{output}\""

rule raxml_check:
    input:
        f"out/{CLUSTER}.msa.fasta"
    output:
        f"out/{CLUSTER}.raxml.rba"
    params:
        cluster_name=CLUSTER
    shell:
        "{RAXML_PATH}/bin/raxml-ng --parse --msa \"{input}\" --model LG+G --prefix \"./out/{params.cluster_name}\""

rule raxml:
    input:
        f"out/{CLUSTER}.raxml.rba"
    output:
        f"out/{CLUSTER}.raxml.bestTree"
    params:
        cluster_name=CLUSTER
    shell:
        "{RAXML_PATH}/bin/raxml-ng --msa \"{input}\" LG+G --prefix \"./out/{params.cluster_name}\""    

rule raxml_stats:
    input:
        f"out/{CLUSTER}.raxml.bestTree"
    shell:
        "python3 scripts/raxml_stats.py \"{input}\""

rule taxtree:
    input:
        f"out/{CLUSTER}.fasta"
    output:
        f"out/{CLUSTER}.tax.tree"
    shell:
        "python3 scripts/taxtree.py \"{input}\""