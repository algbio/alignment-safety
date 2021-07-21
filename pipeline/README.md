## Dependencies
- Diamond   https://github.com/bbuchfink/diamond
- Raxml-ng  https://github.com/amkozlov/raxml-ng
- Muscle    http://www.drive5.com/muscle/
- ete3      https://github.com/etetoolkit/ete
- Hmmer     https://github.com/EddyRivasLab/hmmer

Compile Diamond from source:
```
git clone https://github.com/bbuchfink/diamond
cd diamond-2.0.9
mkdir bin
cd bin
cmake -DEXTRA=ON ..
make
```

Compile Raxml-ng from source:
```
git clone https://github.com/amkozlov/raxml-ng
cd raxml-ng
mkdir build
cd build
cmake ..
make
```

Compile Hmmer and easel from source:
```
git clone https://github.com/EddyRivasLab/hmmer
cd hmmer-3.3.2
./configure --prefix=./
make
cd easel
make
```

ete3 and snakemake via Conda environment
```
conda env create -f environment.yaml
conda activate pipeline
```

---

## How to run the pipeline?

#### 1. Edit `parameters.yaml`:
- `db_file`<br/>
    Path to the protein sequence database (fasta format), e.g. `uniprot_sprot.fasta`.
- `family_db`<br/>
    Path to the Pfam protein domain database: `Pfam-A.seed`.
- `diamond`<br/>
    Path to DIAMOND, e.g. `./../diamond-2.08`.
- `raxml`<br/>
    Path to RaxML-ng, e.g.  `./../raxml-ng`.
- `muscle`<br/>
    Path to Muscle - multiple sequence alignment tool, e.g. `./../muscle`.
- `hmmer`<br/>
    Path to the hmmer programs, e.g. `./../hmmer-3.3.2`.
- `safety`<br/>
    Parent folder of the compiled safety-window program, `./../safety-windows`
- `min_identity`<br/>
    Identity threshhold-% minimum to report and alignment. Used in db clustering process. Tested 20%-90%.
- `sensitivity`<br/>
    `auto`, or any of the DIAMOND sensitivity options: <https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#sensitivity-modes>
- `clustering_algorithm`<br/>
    `mcl` or `multi-step`
- `clustering_min_size`<br/>
    Treshhold of the minimum cluster size to include.<br/>
    Warning: `< 20`, will produce large amount of files
- `clustering_max_size`<br/>
    Treshhold of the maximum cluster size to include.
- `cluster_number`<br/>
    If less than available clusters, `cluster_number` amount of clusters will be chosen randomly to include.<br/>
    Speeds up debugging/testing.
#### 2. Run `separate_clusters` snakemake rule:
- `snakemake -j n separate_clusters`
    `n`: number of parallel processes.
    Clusters the database and separates clusters satisfying the treshholds to `WORK_DIR/`
- `WORK_DIR/fasta/`
    Each fasta-file corresponds to one cluster. Name of the file as well as the first sequence in the file is the reference sequence of that cluster.
- `WORK_DIR/clean/`
    Same as `fasta/`, but fasta sequences are cleaned with no additional information, such as taxonomic id. Needed for some programs, such as Muscle.
- `WORK_DIR/refs/`
    One file containing each cluster's reference sequence in fasta-format. Needed for `phmmer`.
#### 3. Run all rules or some specific rule:
    snakemake -j n rule
- `all`<br/>
    Runs rules: `safe`, `identity`, `hmmsearch`, `hmmscan`, `phmmer` and their prerequisites.
- `safe`<br/>
    Runs Safety-Window-program on all clusters. Outputs to `WORK_DIR/safety/`.
- `identity`<br/>
    Runs `esl-alipid`-program on all clusters to calculate pairwise identities. Outputs to `WORK_DIR/id/`.
- `hmmsearch`<br/>
    Runs `hmmsearch` on all clusters. Outputs to `WORK_DIR/hmmsearch/`.
- `hmmscan`<br/>
    Runs `hmmscan` on all clusters. Outputs to `WORK_DIR/hmmscan/`.
- `phmmer`<br/>
    Runs `phmmer` on all clusters. Outputs to `WORK_DIR/phmmer/`.
    
---

Hmmer documentation: <http://eddylab.org/software/hmmer/Userguide.pdf>
