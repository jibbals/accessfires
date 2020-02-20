#!/bin/bash

## Making videos from figures
#
##Run this from the subfolder with all your images in it


module load ffmpeg
# 
# original ffmpeg call
#ffmpeg -pattern_type glob -framerate 2 -i '*.png' -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p video.mp4


## use optional input arguments for framerate
## ${1:-2} is 2 unless an argument is sent to this script
#ffmpeg -pattern_type glob -framerate ${1:-2} -i '*.png' -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p video.mp4

## AFTER MOVE TO GADI, libx264 encoder no longer works (including profile and crf optional inputs) 
ffmpeg -pattern_type glob -framerate ${1:-2} -i '*.png' -pix_fmt yuv420p video.mp4
