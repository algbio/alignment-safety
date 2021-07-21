## Dependencies
- Diamond   https://github.com/bbuchfink/diamond
- Raxml-ng  https://github.com/amkozlov/raxml-ng
- Muscle    http://www.drive5.com/muscle/
- ete3      https://github.com/etetoolkit/ete

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
