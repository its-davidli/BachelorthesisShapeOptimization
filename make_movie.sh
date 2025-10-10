# Create movie with 
# sudo sh make_movie.sh <folder_name> 

ffmpeg -fflags +genpts -r 1 -i Results/$1/Figures/mesh-%03d.png    -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2"       -vcodec libx264 -y Results/$1/Figures/mesh.mp4
# ffmpeg -r 1 -i TestEField/Results/$1/Figures/state-%03d.png           -vcodec libx264 -y Results/$1/Figures/mesh.mp4
