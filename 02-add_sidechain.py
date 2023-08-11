import mdtraj as md 
import numpy as np
from functions import stack_mdtraj
import os
import sys


seq = sys.argv[1]
straight = bool(int(sys.argv[2]))
anti_parallel = bool(int(sys.argv[3]))
scwrl_path = '/home/gridsan/congwang/scwrl4/Scwrl4'

folder_path_list = [['original', 'straight'], ['parallel', 'anti_parallel']]
folder_name = folder_path_list[0][int(straight)] + '_' + folder_path_list[1][int(anti_parallel)]
output_dir = f'AA/{folder_name}/{seq}'
os.makedirs(output_dir, exist_ok=True)
layer_dist = 0.49


layers_list = []
if straight:
    os.system(f'{scwrl_path} -i backbone/{folder_name}/{seq}/layer_0.pdb -o {output_dir}/layer_0.pdb -h')
    os.system(f'{scwrl_path} -i backbone/{folder_name}/{seq}/layer_1.pdb -o {output_dir}/layer_1.pdb -h')
    for i in range(10):
        if anti_parallel and i % 2 == 1:
            layer_i = md.load(f'{output_dir}/layer_1.pdb')
            layer_i.xyz += np.array([0, 0, layer_dist * (i - 1)])
        else:
            layer_i = md.load(f'{output_dir}/layer_0.pdb')
            layer_i.xyz += np.array([0, 0, layer_dist * i])
        layer_i.save_pdb(f'{output_dir}/layer_{i}.pdb')
        replace_first_line(f'{output_dir}/layer_{i}.pdb')
        layers_list.append(layer_i)
else:
    for i in range(10):
        os.system(f'{scwrl_path} -i backbone/{folder_name}/{seq}/layer_{i}.pdb -o {output_dir}/layer_{i}.pdb -h')
        layer_i = md.load(f'{output_dir}/layer_{i}.pdb')
        replace_first_line(f'{output_dir}/layer_{i}.pdb')
        layers_list.append(layer_i)

AA_all = stack_mdtraj(layers_list)
AA_all.save_pdb(f'{output_dir}/all.pdb')
replace_first_line(f'{output_dir}/all.pdb')

    