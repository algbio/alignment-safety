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
mkdir build
cd build
cmake ..
make
```

Install ete3:
```
pip3 install ete3
```

Before first run:
- Add paths to Snakefile
- After firs run you can comment out line "ncbi.update_taxonomy_database()" in scripts/taxtree.py