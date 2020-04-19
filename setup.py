import subprocess

## Please replace /home/user/ssrbm by your folder path here ##

NAME_FOLDER = '\'/home/user/ssrbm\''

subprocess.call('sed -i "2s@^@rootf=' + NAME_FOLDER + ';@" Align_utils/balign_peps_seed.m', shell=True)
subprocess.call('sed -i "1s@^@rootf=' + NAME_FOLDER + '@" ssRBM.py Align_utils/align_seqpy.py Align_utils/align_to_seedpy.py', shell=True)
