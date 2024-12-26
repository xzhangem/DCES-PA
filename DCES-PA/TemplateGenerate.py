import os
import argparse
import numpy as np
import utils.input_output as io
import plotly.graph_objects as go
from H2_match import H2StandardIterative, H2MultiRes
from varifold_mean_pick import varifold_mean_file
import torch
import trimesh
import H2_stats as stats


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Template Generation')
    parser.add_argument('--data_directory', type=str,
            help='file name of mesh data')
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
    parser.add_argument('--temp_save', type=str, default='./save_template_file/',
            help='file saves the mesh template')

    args = parser.parse_args()
    if not os.path.exists(args.temp_save):
        os.makedirs(args.temp_save)

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

    initial_template = varifold_mean_file(mesh_file=args.data_directory)
    initial_template = trimesh.load(args.data_directory + initial_template)
    initial_template_np = [np.array(initial_template.vertices)/5, np.array(initial_template.faces)]

    mesh_samples = []
    mesh_dir = os.listdir(args.data_directory)
    for i in range(len(mesh_dir)):
        print("\n Sample {}/{} pre-align with template".format(i+1, len(mesh_dir)))
        mesh = trimesh.load(args.data_directory + mesh_dir[i])
        trans_matrix, cost = trimesh.registration.mesh_other(mesh, initial_template, samples=2000, scale=True, icp_first=20, icp_final=60)
        trans_v = trimesh.transformations.transform_points(mesh.vertices, trans_matrix)
        trans_v = trans_v / 5
        mesh_samples.append([np.array(trans_v), np.array(mesh.faces)])

    V_mu, F_mu = stats.H2UnparamKMean(mesh_samples, initial_template_np, args.a0, args.a1, args.b1, args.c1, args.d1, args.a2, 
            paramlist, N=None, geods_len=args.geods_len)
    io.saveData(file_name=args.temp_save+'karcher_mean', extension='ply', V=V_mu*5, F=F_mu, Rho=None, color=None)
