#!/bin/bash
# to generate Dockerfile
# pip install neurodocker # only the first time, then:
# neurodocker generate docker --base-image debian:buster \
# --pkg-manager apt \
# --ants version=2.4.1 \
# --fsl version=6.0.5 > Dockerfile
# after that lines were commented out and ants version was replaced by 2.4.3 throughout the resulting Dockerfile

# run these the first time:
# docker buildx create --name mstools
# docker builx use mstools
docker buildx build --platform linux/amd64 \
    --build-arg BUILDKIT_INLINE_CACHE=1 \
    --cache-from=pennsive/mstools:4.2 \
    --cache-to=pennsive/mstools:4.2 \
    -t pennsive/mstools:4.2 \
    --push . &> $(date +"%Y-%m-%d_%T")_build.log
       