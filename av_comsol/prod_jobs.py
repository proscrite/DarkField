import os
from time import sleep
exec_template_file = '/scratch/herrero/el/jobs/job_template.sh'
exec_template = open(exec_template_file).read()

JOBSDIR = '.'
exec_params = {'jobsdir' : JOBSDIR}

for j in range(10, 110, 10):
    jobfilename = JOBSDIR + '/' + 'el_ar_{}.sh'.format(j)
    jobfile = open(jobfilename, 'w')
    jobfile.write(exec_template.format(**exec_params))
    cmd = './el {} 300 > Results/ar/{}/el_out_300.dat'.format(j,j)
    jobfile.write(cmd)
jobfile.close()

#send jobs
for i in range(10, 110, 10):
    cmd = 'qsub {}/el_ar_{}.sh'.format(JOBSDIR, i)
    print(cmd)
    os.system(cmd)
    sleep(0.5)
