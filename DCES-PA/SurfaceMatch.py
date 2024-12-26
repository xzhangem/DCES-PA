import os
import argparse
import numpy as np
import utils.input_output as io
import plotly.graph_objects as go
from H2_match import H2StandardIterative, H2MultiRes
import torch 
import trimesh
from scipy.io import savemat

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Matching template to the current sample')
    parser.add_argument('--template_name', type=str,
            help='template name, end with .ply')
    parser.add_argument('--data_directory', type=str,
            help='file name of mesh data')
    parser.add_argument('--pre_align', type=bool, default=False,
            help='pre-align template with sample')
    parser.add_argument('--a0', type=float, default=1,
            help='a0 coeffecient')
    parser.add_argument('--a1', type=float, default=200,
            help='a1 coeffecient')
    parser.add_argument('--b1', type=float, default=0,
            help='b1 coeffecient')
    parser.add_argument('--c1', type=float, default=200,
            help='c1 coeffecient')
    parser.add_argument('--d1', type=float, default=0,
            help='d1 coeffecient')
    parser.add_argument('--a2', type=float, default=200,
            help='a2 coeffecient')
    parser.add_argument('--geods_len', type=int, default=2,
            help='numbers of state for approximating geodesics')
    parser.add_argument('--post_align', type=bool, default=False,
            help='post-align result with sample')
    parser.add_argument('--mesh_savefile', type=str, default='./save_mesh_file/',
            help='file saves the mesh result')
    parser.add_argument('--mat_savefile', type=str, default='./save_mat_file/',
            help='file saves the mat result')

    args = parser.parse_args()
    if not os.path.exists(args.mesh_savefile):
        os.makedirs(args.mesh_savefile)
    if not os.path.exists(args.mat_savefile):
        os.makedirs(args.mat_savefile)

    param1_mul = { 'weight_coef_dist_T': 10**1,'weight_coef_dist_S': 10*2.0,'sig_geom': 5,\
            'max_iter': 5000,'time_steps': 2, 'tri_unsample': True, 'index':0}
    param2_mul = {'weight_coef_dist_T': 10**1,'weight_coef_dist_S': 10*2.0,'sig_geom': 2.5,\
            'max_iter': 5000,'time_steps': 2, 'tri_unsample': False, 'index':1}
    param3_mul =  {'weight_coef_dist_T': 10**1,'weight_coef_dist_S': 10*2.0,'sig_geom':1.0,\
            'max_iter': 5000, 'time_steps': 2, 'tri_unsample': True, 'index':1}
    param4_mul = {'weight_coef_dist_T': 10**1,'weight_coef_dist_S': 10*2.0,'sig_geom': 0.5,\
            'max_iter': 5000, 'time_steps': 2, 'tri_unsample': False, 'index':2}
    param5_mul = { 'weight_coef_dist_T': 10**1,'weight_coef_dist_S': 10*2.0,'sig_geom':0.1,\
            'max_iter': 5000, 'time_steps': 2, 'tri_unsample': False, 'index':2}

    paramlist = [param1_mul, param2_mul, param3_mul, param4_mul, param5_mul]

    template_mesh = trimesh.load(args.template_name)

    mesh_samples = []
    mesh_dir = os.listdir(args.data_directory)
    for i in range(len(mesh_dir)):
        mesh = trimesh.load(args.data_directory + mesh_dir[i])
        mesh_samples.append(mesh)

    for i in range(len(mesh_dir)):
        print("Processing sample of {}/{}".format(i+1, len(mesh_dir)))
        mesh = mesh_samples[i]
        if args.pre_align == True:
            trans_matrix, cost = trimesh.registration.mesh_other(template_mesh, mesh, samples=2000, scale=True, icp_first=20, icp_final=60)
            trans_v = trimesh.transformations.transform_points(template_mesh.vertices, trans_matrix)
            trans_v = trans_v / 10
        else:
            trans_v = template_mesh.vertices / 10
        template_array = [np.array(trans_v), np.array(template_mesh.faces)]
        sample_array = [np.array(mesh.vertices)/10, np.array(mesh.faces)]

        geod, F0 = H2MultiRes(template_array, sample_array, args.a0, args.a1, args.b1, args.c1, args.d1, args.a2, 2, paramlist)

        if args.post_align == True:
            result_mesh = trimesh.Trimesh(vertices=geod[-1]*10, faces=F0)
            trans_matrix, cost = trimesh.registration.mesh_other(result_mesh, mesh, samples=2000, scale=True, icp_first=20, icp_final=60)
            trans_v = trimesh.transformations.transform_points(result_mesh.vertices, trans_matrix)
        else:
            trans_v = vertices=geod[-1]*10

        io.saveData(file_name= args.mesh_savefile+mesh_dir[i].split('.')[0], extension='ply', V=trans_v, F=F0, Rho=None, color=None)
        savemat(args.mat_savefile+mesh_dir[i].split('.')[0]+'.mat', {'V':trans_v, 'F':F0+1})

