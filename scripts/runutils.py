import os


def runcmd(folder, cmd):
    command = ""
    fname = ""
    for arg in cmd:
        if ' ' in arg:
            command += ("\"" + arg + "\"")
            fname += arg.replace(' ', '_')
        else:
            command += arg
            fname += arg
        command += " "
        fname += "_"
    command += "> "
    fname = fname[:-1]
    fname = folder + "/" + fname
    command += fname
    command = "./" + command
    print(command)
    os.system(command)
    analys = "python3 ./runtest.py cat " + fname + " > " + fname + ".stat"
    print(analys)
    os.system(analys)
