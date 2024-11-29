import sys, os

#doing trigCorex instead of plot_2mu2e?
trig = False

#use backup file instead of normal one?
backup = False

#the errfiles.txt should be the result of this command:
# grep "OSError" logs/*.err >errfiles.txt
inname = "errfiles.txt" #"missing_files.txt"
if backup:
    inname = "errfiles_backup.txt"

infile = open(inname, "r")

old_let = ""
old_setn = -1
old_num = -1
for line in infile:
    words = line.strip().split("/")
    parts = words[1].split("_")

    direc = parts[1]
    if "test" in direc:
        direc = parts[2]
        num = int(parts[3].split('.')[0])
    else:
        num = int(parts[2].split('.')[0])
    let = direc[0]
    setnum = direc[1]
    if setnum == 'v':
        let = let + setnum + direc[2]
        setnum = direc[3]
    setn = int(setnum)
    
    #don't submit the job a 2nd time if it had 2 lines of errors.
    if let == old_let and setn == old_setn and num == old_num:
        continue

    if trig:
        template = "sub_trigCorex_template"
    else:
        template = "sub_template"
    template += "%s.jdl"%("_backup" if backup else "") 
    #get the new jdl file
    #if backup:
    os.system("sed \"s/letter=C/letter=%s/g\" %s >temp0.jdl"%(let, template))
    #else:
    #    os.system("sed \"s/letter=C/letter=%s/g\" sub_template.jdl >temp0.jdl"%(let))
    os.system("sed \"s/setnum=0/setnum=%d/g\" temp0.jdl >temp1.jdl"%(setn))
    os.system("sed \"s/\$(PROCESS)/%d/g\" temp1.jdl >temp2.jdl"%(num))
    if backup:
        os.system("sed \"s/Queue 32/Queue 1/g\" temp2.jdl > sub_jobs_backup.jdl") 
        os.system("condor_submit sub_jobs_backup.jdl")
    else:
        os.system("sed \"s/Queue 32/Queue 1/g\" temp2.jdl > sub_jobs.jdl") 
        os.system("condor_submit sub_jobs.jdl")

    old_let = let
    old_setn = setn
    old_num = num


