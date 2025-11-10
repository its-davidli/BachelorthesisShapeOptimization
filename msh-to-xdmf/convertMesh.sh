# Convert  .geo into .xdmf mesh files 
# sh msh-to-xdmf/convertMesh.sh <meshfile_name> <subfolder_name> <dimension>

gmsh -$3 msh-to-xdmf/$1.geo
python3 msh-to-xdmf/convert_msh_to_xdmf.py $1 $2 $3 
