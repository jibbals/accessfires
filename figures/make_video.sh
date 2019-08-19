#!/bin/bash

# Run this from the subfolder with all your images in it
module load ffmpeg

ffmpeg -pattern_type glob -framerate 2 -i 'fig_201*.png' -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p video.mp4
