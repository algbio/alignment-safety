#!/usr/bin/python3

"""
Running simple calcultion experiments for emerald.
Outputting a csv file things such as coverage and number
of safety windows, and time and memory consumption.
Creates some log files in ./experiments directories as well.
Assumes the program ./main is in the running directory.
"""

import sys
import os
import platform

# Helper command for run
def time_cmd(cmd, log):
    if platform.system() == "Darwin":
        # Mac OS X System. Install GNU time with homebrew first
        return f"/opt/homebrew/opt/gnu-time/libexec/gnubin/time -o {log} -a -v {cmd}"
    # GNU Linux
    return f"/usr/bin/time -o {log} -a -v {cmd}"

# runner function
def run(cmd, log="", time=True, verbose=True):
    if verbose:
        print(cmd)
    if log != "":
        open(log, 'a').write(cmd + '\n')
        if time:
            cmd = time_cmd(cmd, log)
    try:
        os.system(cmd)
    except KeyboardInterrupt:
        exit(1)

# Startup check
if len(sys.argv) != 3 or sys.argv[1] != "-p":
    print("How to run: python3", sys.argv[0], "-p {fasta_file}")
    exit(1)

# path contains all the gfa files
path = sys.argv[2]
direc = os.fsencode(path)

exp_path = os.path.abspath('.') + "/experiments/"
if not os.path.exists('experiments'):
    os.makedirs('experiments')

# Writing csv file
sp = path.rsplit('/')[-1]
csv = exp_path + f"{sp}.csv"
print(f"Writing to csv file {csv}")
with open(csv, 'w') as f: # use mode 'w' to overwrite last file
    f.write("Alpha,Delta,AVG number of safety windows per sequence,STD Dev of number of safety windows per sequence,Coverage,Runtime,Memory")

# Content written into csv file will first be saved into this list
csvlines = []

def time_to_csv(logfile):
    # TODO
    return

def output_to_csv(sw_out):
    # TODO
    return

# Run with all (Alpha, Delta) inputs
inputs = [(0.5, 0), (0.7, 0), (0.99, 0), (0.5, 2), (0.7, 2), (0.99, 2), (0.5, 5), (0.7, 5)]
prg = f"./main -f {path}"
logfile = f"{exp_path}/calc.log"
sw_out = f"{exp_path}/{sp}.out"
for inp in inputs:
    run(prg + f" -a {inp[0]} -d {inp[1]} > {sw_out}", logfile)
    time_to_csv(logfile)
    output_to_csv(sw_out)


with open(csv, 'a') as f:
    for line in csvlines:
        for i in range(len(line)):
            f.write(f"{line[i]}")
            f.write("\n") if i + 1 == len(line) else f.write(",")
