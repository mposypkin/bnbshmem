import os
import subprocess
import sys
import csv
import statistics

def mmdev(serl):
     m = statistics.mean(serl)
     minv = min(serl)
     maxv = max(serl)
     return max(m - minv, maxv - m)

def makeStat(name, serl):
    print(name, " mean: ", statistics.mean(serl))
    print(name, " stddev: ", statistics.stdev(serl))
    print(name, " min: ", min(serl))
    print(name, " max: ", max(serl))
    print(name, " maxmindev: ", mmdev(serl))


output = subprocess.run(sys.argv[1:],stdout=subprocess.PIPE)
sout = output.stdout.decode()
p = sout.find('Statistics:')
ls = sout[p:].splitlines()
cls = ls[1:]

rd = csv.reader(cls)
timel = []
stepl = []
tpsl = []
for row in rd:
    timel.append(float(row[0]))
    stepl.append(float(row[1]))
    tpsl.append(float(row[2]))
    #print(row)

makeStat("Time", timel)
makeStat("Steps", stepl)
makeStat("TPS", tpsl)
#print(output.stdout.decode())
#process = os.popen(sys.argv[1:])
#strout = process.read()
#process.close()
#print(strout)
