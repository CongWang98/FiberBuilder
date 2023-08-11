import mdtraj as md 
from functions import rotate_molecule, norm_vector, stack_mdtraj
import numpy as np
import os
import sys


seq = sys.argv[1]
if seq == 'WRWW':
    res_list = ['TRP', 'ARG', 'TRP', 'TRP']
elif seq == 'WWWW':
    res_list = ['TRP', 'TRP', 'TRP', 'TRP']
res_idx_list = [2, 3, 4, 5]
res_name_mapping = dict(zip(res_idx_list, res_list))


straight = bool(int(sys.argv[2]))
anti_parallel = bool(int(sys.argv[3]))
folder_path_list = [['original', 'straight'], ['parallel', 'anti_parallel']]
folder_name = folder_path_list[0][int(straight)] + '_' + folder_path_list[1][int(anti_parallel)]
output_dir = f'backbone/{folder_name}/{seq}'
os.makedirs(output_dir, exist_ok=True)


fiber = md.load('start_fiber.pdb')
backbone_idx = fiber.top.select('backbone')
fiber_backbone = fiber.atom_slice(backbone_idx)


table, bonds = fiber_backbone.topology.to_dataframe()
for res_idx, new_residue_name in res_name_mapping.items():
    table.loc[table['resSeq'] == res_idx, 'resName'] = new_residue_name
new_topology = md.Topology.from_dataframe(table, bonds)
fiber_backbone.topology = new_topology


layers_list = []
for i in range(10):
    layer_i_idx = fiber_backbone.top.select(f'chainid {i*7} to {i*7 + 6}')
    layer_i = fiber_backbone.atom_slice(layer_i_idx)
    layer_i_center = layer_i.xyz[0].mean(axis=0)
    layers_list.append((layer_i, layer_i_center))
layers_list.sort(key=lambda x: x[1][2])


layer_dis = []
for i in range(len(layers_list) - 1):
    layer_dis.append(np.linalg.norm(layers_list[i][1] - layers_list[i+1][1]))
print(layer_dis)
print(np.mean(layer_dis))
layer_dist = 0.49


layers_list = [layer_i[0] for layer_i in layers_list]
layers_list[0].save_pdb(f'{output_dir}/_ref.pdb')

for i, layer_i in enumerate(layers_list):
    if straight:
        layer_i = md.load(f'{output_dir}/_ref.pdb')
        layer_i.xyz += np.array([0, 0, layer_dist * i])
    if anti_parallel and i % 2 == 1:
        rotate_list = []
        layer_i.save_pdb(f'{output_dir}/layer_{i}_before_rotate.pdb')
        for j in range(7):
            chain_ij_idx = layer_i.top.select(f'chainid {j}')
            chain_ij = layer_i.atom_slice(chain_ij_idx)
            chain_ij_mean = chain_ij.xyz[0].mean(axis=0, keepdims=True)
            chain_ij.xyz[0] -= chain_ij_mean
            point0 = np.array([0, 0, 0])
            point1 = chain_ij.xyz[0][0]
            point2 = chain_ij.xyz[0][-1]
            rotate_axis = norm_vector(point0, point1, point2)
            chain_ij_rotated = rotate_molecule(chain_ij, rotate_axis, np.pi)
            chain_ij_rotated.xyz[0] += chain_ij_mean
            rotate_list.append(chain_ij_rotated)
        layer_i = stack_mdtraj(rotate_list)
    layer_i.save_pdb(f'{output_dir}/layer_{i}.pdb')
    layers_list[i] = layer_i

fiber_backbone_all = stack_mdtraj(layers_list)
fiber_backbone_all.save_pdb(f'{output_dir}/all.pdb')

    