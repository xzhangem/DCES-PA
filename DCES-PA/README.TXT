The implementation of "DCES-PA: Deformation-controllable elastic shape model for 3D bone proliferation analysis using hand HR-pQCT images"

The Python enviornment can be referred to enviornment.yml 

To generate template joint mesh:

	python TemplateGenerate.py --data_directory '../MCP_mesh_volume/MCP2_health_mseh/'

To match the generated template to the given sample for proliferation analysis:

	
	python SurfaceMatch.py --template_name 'save_template_file/karcher_mean.ply' --data_directory '../MCP_mesh_volume/MCP2_PsA_mesh/'

For visualization:

	DCES_display.m


Ackknowledgement: The shape analysis framework is built and modified based on "Hartman et al., Elastic shape analysis of surfaces with second-order Sobolev metrics: a comprehensive numerical framework. IJCV 2023".
The extrinsically-weighted second-order metric is included in 'enr/DDG.py'.  