#!/bin/bash

DIR=$1

if [[ "$DIR" == "" ]]
then
  DIR=sapling
fi

ffmpeg -r 30 -s 1280x720  -i image.%05d.tga -vcodec libx264 -crf 15 vid.mp4

