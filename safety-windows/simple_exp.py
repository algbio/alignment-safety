#!/usr/bin/python3

"""
Running simple calcultion experiments for emerald.
Outputting a csv file things such as coverage and number
of safety windows, and time and memory consumption.
Creates some log files in ./experiments directories as well.
Assumes the program ./main is in the running directory.
"""

import math
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
    f.write("Alpha,Delta,AVG number of safety windows per sequence,STD Dev of number of safety windows per sequence,Coverage,Runtime,Memory\n")

# Content written into csv file will first be saved into this list
csvlines = []

def time_to_csv(logfile):
    with open(logfile, 'r') as lgf:
        lines = lgf.readlines()[-23:]

        runtime = float(lines[1].split(' ')[-1]) + float(lines[2].split(' ')[-1])
        mem = int(lines[9].split(' ')[-1])
        csvlines[-1].extend([runtime, mem])

def output_to_csv(sw_out):
    with open(sw_out, 'r') as f:
        lines = f.readlines()
        ref = lines[1]
        n = len(ref)
        cov = [False] * n
        COV = 0
        AVG = 0
        STDDEV = 0
        swss = []

        am = int(lines[2])
        j = 3
        for i in range(am - 1):
            sws = int(lines[j + 2])
            swss.append(sws)
            AVG += sws
            for k in range(sws):
                L, R, foo1, foo2 = map(int, lines[j + k + 3].split())
                for i in range(L, R + 1):
                    cov[i] = True
            j += sws + 3
        AVG *= 1.0/am

        for b in cov:
            if b:
                COV += 1
        COV *= 1.0/n

        for sws in swss:
            STDDEV += (sws - AVG) * (sws - AVG)
        STDDEV = math.sqrt(STDDEV)
        if am > 1:
            STDDEV /= n - 1

        csvlines[-1].extend([AVG, STDDEV, COV])



# Run with all (Alpha, Delta) inputs
inputs = [(0.5, 0), (0.7, 0), (0.99, 0), (0.5, 2), (0.7, 2), (0.99, 2), (0.5, 5), (0.7, 5)]
prg = f"./main -f {path}"
logfile = f"{exp_path}/calc.log"
sw_out = f"{exp_path}/{sp}.out"
for inp in inputs:
    run(prg + f" -a {inp[0]} -d {inp[1]} > {sw_out}", logfile)
    csvlines.append([inp[0], inp[1]])
    output_to_csv(sw_out)
    time_to_csv(logfile)


with open(csv, 'a') as f:
    for line in csvlines:
        for i in range(len(line)):
            f.write(f"{line[i]}")
            f.write("\n") if i + 1 == len(line) else f.write(",")
