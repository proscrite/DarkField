import os
from time import sleep
exec_template_file = '/scratch/herrero/avComsol/jobs/atlas_job_template.sh'
exec_template = open(exec_template_file).read()

JOBSDIR = '.'
exec_params = {'jobsdir' : JOBSDIR}

sym_wire = ['62dot5', '75', '87dot5', '100', '112dot5', '125']

#write .sh jobs
for j in sym_wire:
    jobfilename = JOBSDIR + '/' + 'av_sym_{}um.sh'.format(j)
    jobfile = open(jobfilename, 'w')
    jobfile.write(exec_template.format(**exec_params))
    cmd = './el {}um 500 > results/sym/{}/av_sym_500_out.dat'.format(j,j)
    jobfile.write(cmd)
jobfile.close()

#send jobs
for i in sym_wire:
    cmd = 'qsub {}/av_sym_{}um.sh'.format(JOBSDIR, i)
    print(cmd)
    os.system(cmd)
    sleep(0.5)
