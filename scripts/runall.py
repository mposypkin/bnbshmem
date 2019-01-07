import os
import time
import runutils

def run(folder, bname, nruns, rec, eps, nsteps):
    runutils.runcmd(folder, ["bnbomp.exe", nruns, bname, rec, eps, nsteps, "12"])
    time.sleep(30)
    runutils.runcmd(folder, ["gpamigo.exe", nruns, bname, rec, eps, nsteps, "12"])
    time.sleep(30)
    runutils.runcmd(folder, ["lpamigo.exe", nruns, bname, rec, eps, nsteps, "12"])
    time.sleep(30)
    runutils.runcmd(folder, ["bnbatomic.exe", nruns, bname, rec, eps, nsteps, "4096", "1000"])

dir = "./results"
if not os.path.exists(dir):
    os.makedirs(dir)

# run(dir, "Cluster2D2 function", "32", "uknrec", "0.1", "100000000")
# run(dir, "Hartman 6 function", "32", "uknrec", "0.01", "100000000")
# run(dir, "Biggs EXP6 Function", "32", "uknrec", "0.1", "10000000")
runutils.runcmd(dir, ["bnbatomic.exe", "32", "Biggs EXP6 Function", "uknrec", "0.1", "10000000", "4096", "1000"])
