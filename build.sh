#!/bin/bash
# to generate Dockerfile
# pip install neurodocker
# neurodocker generate docker --base-image debian:buster \
# --pkg-manager apt \
# --ants version=2.4.1 \
# --fsl version=6.0.5 > Dockerfile
# then lines were commented out and ants version was replaced by 2.4.3 throughout
docker build -t pennsive/mstools:base .