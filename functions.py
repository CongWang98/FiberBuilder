import numpy as np
import mdtraj as md



def norm_vector(point1, point2, point3):
    normal_vector = np.cross(point2 - point1, point3 - point1)
    normal_vector /= np.linalg.norm(normal_vector)
    return normal_vector


def rotation_matrix(axis, theta):
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    cross_matrix = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])
    rotation_matrix = (
        cos_theta * np.identity(3) +
        sin_theta * cross_matrix +
        (1 - cos_theta) * np.outer(axis, axis)
    )
    return rotation_matrix


#def norm_vector_atom(traj, atom_idx1, atom_idx2, atom_idx3):
#    positions = traj.xyz[0][[atom_idx1, atom_idx2, atom_idx3]]
#    point1, point2, point3 = positions[0], positions[1], positions[2]
#    vector1 = point2 - point1
#    vector2 = point3 - point1
#    normal_vector = np.cross(vector1, vector2)
#    return normal_vector


def rotate_molecule(traj, axis, theta):
    rot_matrix = rotation_matrix(axis, theta)
    rotated_trajectory = []
    for frame in traj.xyz:
        rotated_frame = np.dot(frame, rot_matrix.T)
        rotated_trajectory.append(rotated_frame)
    rotated_trajectory = np.array(rotated_trajectory)
    rotated_traj = md.Trajectory(rotated_trajectory, topology=traj.topology)
    return rotated_traj


def stack_mdtraj(traj_list):
    merged_traj = traj_list[0]
    for i in range(1, len(traj_list)):
        merged_traj = merged_traj.stack(traj_list[i])
    return merged_traj


def replace_first_line(file_name):
    '''
    Mdtraj saves the first line of the pdb file as: MODEL        0
    Martinize2 has a stupid bug that it cannot recognize this line and its developers havent fixed it for 3 years!
    This function replaces the first line MODEL        0 with MODEL        1 so that martinize2 can read it.
    '''
    with open(file_name, 'r') as f:
        lines = f.readlines()
    lines[0] = 'MODEL        1\n'
    with open(file_name, 'w') as f:
        f.writelines(lines)