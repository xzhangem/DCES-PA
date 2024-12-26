clear all;
close all;
clc

bone_file = ''; % The .mat format file of the joint bone mesh 
syn_file = ''; % The .mat format file of the matching result 

bone_dir = dir(fullfile(syn_file, '*.mat'));
bone_num = length(bone_dir);
th = 3.5; % above threshold part may be taken as proliferations
db_radium = 10; % parameters for DBSCAN
point_th = 35; % Proliferations formed by more than point_th points may be possible ones

for i = 1 : bone_num
    bone_name = [bone_file, bone_dir(i).name];
    syn_name = [syn_file, bone_dir(i).name];
    load(syn_name)
    faces = int32(F);
    vertices = V;
    FV.faces = faces;
    FV.vertices = vertices;
    
    load(bone_name)
    points = V;
    [distances, surface_points] = point2trimesh(FV, 'QueryPoints', points, 'UseSubSurface', false);
    
    point_num = length(distances);
    nonzero_distances = single(distances >= th);
    nonzero_points = points .* nonzero_distances;
    nonzero_match_points = surface_points .* nonzero_distances;
    
    nonzero_points(any(nonzero_points, 2)==0, :) = [];
    nonzero_match_points(any(nonzero_match_points, 2)==0, :) = [];
    diseased_points = [nonzero_points; nonzero_match_points];
    
    % Clustering possible proliferation regions with DBSCAN
    idx = dbscan(diseased_points, db_radium, 1);
    idx_unique = unique(idx);
    if min(idx_unique) == -1
        idx_unique = idx_unique(2:end, :);
    end
    vol_vec = zeros(length(idx_unique), 1);
    
    % Alpha-shape for proliferation
    K_cell = {};
    ost_points_cell = {};
    
    for k = 1 : length(idx_unique)
        ost_points = diseased_points .* single(idx == k);
        ost_points(any(ost_points, 2)==0,:) = [];
        [K, vol] = boundary(double(ost_points), 0.5);
        vol_vec(k, 1) = vol;
        K_cell{k} = K;
        ost_points_cell{k} = ost_points;
    end
    
    % Visualize the possible proliferations
    for n = 1 : length(vol_vec)
        K = K_cell{n};
        ost_points = ost_points_cell{n};
        d = zeros(size(ost_points(:,1)));
        [ost_point_num, ~] = size(ost_points_cell{n});
        if ost_point_num >= point_th
            hold on
            figure
            trisurf(K, ost_points(:,1), ost_points(:,2), ost_points(:,3), 'FaceAlpha', 0.8);
        end
    end
    
    % Display the part with positive distances indicating proliferations
    distances = single(distances > 0) .* distances;
    distances = distances * 82 / 1000;
    
    hold on
    figure
    patch(FV,'Facecolor','c','FaceAlpha',.5); xlabel('x'); ylabel('y'); 
    zlabel('z'); axis equal; hold on; grid on
    plot3M = @(XYZ,varargin) plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),varargin{:});
    plot3M(points,'*r')
    plot3M(surface_points,'*k')
    plot3M(reshape([shiftdim(points,-1);shiftdim(surface_points,-1);
    shiftdim(points,-1)*NaN],[],3),'k'); 
    grid on
    abs_distances = abs(distances);
    max_abs_distances = max(abs_distances);
    dir_dist=(1).*distances;
    
    figure
    trisurf(F,V(:,1),V(:,2),V(:,3), dir_dist(:,1), 'Edgecolor', 'none') %dir_dist(:,1),
    %trisurf(recon_mesh_F, recon_mesh_V(:,1), recon_mesh_V(:,2), recon_mesh_V(:,3), dir_dist(:,1), 'Edgecolor', 'none')
    axis equal
    material dull
    lightangle(-45,30)
    shading interp
    FaceLighting = 'gouraud';
    AmbientStrength = 0.9;
    DiffuseStrength = 0.8;
    SpecularStrength = 0.9;
    SpecularExponent = 25;
    BackFaceLighting = 'unlit';
    caxis([-max_abs_distances max_abs_distances]);
    cob = colorbar('Ticks',-max_abs_distances:5*82/1000:max_abs_distances);
    cob.Label.String = 'Millimeters';
    cob.Label.FontSize = 18;
    colormap(jet)
    colormapeditor
end
