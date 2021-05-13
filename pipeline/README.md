## Dependencies
- Diamond   https://github.com/bbuchfink/diamond
- Raxml-ng  https://github.com/amkozlov/raxml-ng
- Muscle    http://www.drive5.com/muscle/
- ete3      https://github.com/etetoolkit/ete

Install Diamond:
```
git clone https://github.com/bbuchfink/diamond
cd diamond-2.0.9
mkdir bin
cd bin
cmake -DEXTRA=ON ..
make
```

Install Raxml-ng:
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

Before first run:
- Add paths to Snakefile