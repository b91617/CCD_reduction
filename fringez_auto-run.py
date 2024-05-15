#%%
import subprocess
import os

# Run Fringez
path = f'{os.getcwd()}/LOT/Test'

for folder in os.listdir(path):
    path_list = os.path.join(path, folder)
    if os.path.isdir(path_list) and folder != '._':
        working_dir = os.path.join(path, folder) + '/calibrated_data/reduced_lights'
        print(working_dir)

        command = ['fringez-generate',
                   'fringez-clean --all-images-in-folder --fringe-model-folder="{}"'.format(working_dir)]
    
        for cmd in command:
            proc = subprocess.Popen(cmd,
                                    shell=True,
                                    stdout=subprocess.PIPE,
                                    cwd=working_dir)
    
            for line in iter(proc.stdout.readline, b''):
                print(line.decode('utf-8'), end='')
    
            proc.wait()
            

