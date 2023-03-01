 # Generated by Neurodocker and Reproenv.

FROM debian:buster
ENV ANTSPATH="/opt/ants-2.4.3/" \
    PATH="/opt/ants-2.4.3:$PATH"
RUN apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
           ca-certificates \
           curl \
           unzip \
    && rm -rf /var/lib/apt/lists/* \
    && echo "Downloading ANTs ..." \
    && curl -fsSL -o ants.zip https://github.com/ANTsX/ANTs/releases/download/v2.4.3/ants-2.4.3-centos7-X64-gcc.zip \
    && unzip ants.zip -d /opt \
    && mv /opt/ants-2.4.3/bin/* /opt/ants-2.4.3 \
    && rm ants.zip
ENV FSLDIR="/opt/fsl-6.0.5" \
    PATH="/opt/fsl-6.0.5/bin:$PATH" \
    FSLOUTPUTTYPE="NIFTI_GZ" \
    FSLMULTIFILEQUIT="TRUE" \
    FSLTCLSH="/opt/fsl-6.0.5/bin/fsltclsh" \
    FSLWISH="/opt/fsl-6.0.5/bin/fslwish" \
    FSLLOCKDIR="" \
    FSLMACHINELIST="" \
    FSLREMOTECALL="" \
    FSLGECUDAQ="cuda.q"
RUN apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
           bc \
           ca-certificates \
           curl \
           dc \
           file \
           libfontconfig1 \
           libfreetype6 \
           libgl1-mesa-dev \
           libgl1-mesa-dri \
           libglu1-mesa-dev \
           libgomp1 \
           libice6 \
           libopenblas-base \
           libxcursor1 \
           libxft2 \
           libxinerama1 \
           libxrandr2 \
           libxrender1 \
           libxt6 \
           nano \
           sudo \
           wget \
    && rm -rf /var/lib/apt/lists/* \
    && echo "Downloading FSL ..." \
    && mkdir -p /opt/fsl-6.0.5 \
    && curl -fL https://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-6.0.5-centos7_64.tar.gz \
    | tar -xz -C /opt/fsl-6.0.5 --strip-components 1 
    #&& echo "Installing FSL conda environment ..." \
    #&& bash /opt/fsl-6.0.5/etc/fslconf/fslpython_install.sh -f /opt/fsl-6.0.5

# Save specification to JSON.
RUN printf '{ \
  "pkg_manager": "apt", \
  "existing_users": [ \
    "root" \
  ], \
  "instructions": [ \
    { \
      "name": "from_", \
      "kwds": { \
        "base_image": "debian:buster" \
      } \
    }, \
    { \
      "name": "env", \
      "kwds": { \
        "ANTSPATH": "/opt/ants-2.4.3/", \
        "PATH": "/opt/ants-2.4.3:$PATH" \
      } \
    }, \
    { \
      "name": "run", \
      "kwds": { \
        "command": "apt-get update -qq\\napt-get install -y -q --no-install-recommends \\\\\\n    ca-certificates \\\\\\n    curl \\\\\\n    unzip\\nrm -rf /var/lib/apt/lists/*\\necho \\"Downloading ANTs ...\\"\\ncurl -fsSL -o ants.zip https://github.com/ANTsX/ANTs/releases/download/v2.4.3/ants-2.4.3-centos7-X64-gcc.zip\\nunzip ants.zip -d /opt\\nmv /opt/ants-2.4.3/bin/* /opt/ants-2.4.3\\nrm ants.zip" \
      } \
    }, \
    { \
      "name": "env", \
      "kwds": { \
        "FSLDIR": "/opt/fsl-6.0.5", \
        "PATH": "/opt/fsl-6.0.5/bin:$PATH", \
        "FSLOUTPUTTYPE": "NIFTI_GZ", \
        "FSLMULTIFILEQUIT": "TRUE", \
        "FSLTCLSH": "/opt/fsl-6.0.5/bin/fsltclsh", \
        "FSLWISH": "/opt/fsl-6.0.5/bin/fslwish", \
        "FSLLOCKDIR": "", \
        "FSLMACHINELIST": "", \
        "FSLREMOTECALL": "", \
        "FSLGECUDAQ": "cuda.q" \
      } \
    }, \
    { \
      "name": "run", \
      "kwds": { \
        "command": "apt-get update -qq\\napt-get install -y -q --no-install-recommends \\\\\\n    bc \\\\\\n    ca-certificates \\\\\\n    curl \\\\\\n    dc \\\\\\n    file \\\\\\n    libfontconfig1 \\\\\\n    libfreetype6 \\\\\\n    libgl1-mesa-dev \\\\\\n    libgl1-mesa-dri \\\\\\n    libglu1-mesa-dev \\\\\\n    libgomp1 \\\\\\n    libice6 \\\\\\n    libopenblas-base \\\\\\n    libxcursor1 \\\\\\n    libxft2 \\\\\\n    libxinerama1 \\\\\\n    libxrandr2 \\\\\\n    libxrender1 \\\\\\n    libxt6 \\\\\\n    nano \\\\\\n    sudo \\\\\\n    wget\\nrm -rf /var/lib/apt/lists/*\\necho \\"Downloading FSL ...\\"\\nmkdir -p /opt/fsl-6.0.5\\ncurl -fL https://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-6.0.5-centos7_64.tar.gz \\\\\\n| tar -xz -C /opt/fsl-6.0.5 --strip-components 1" \
      } \
    } \
  ] \
}' > /.reproenv.json
# End saving to specification to JSON.

ARG R_VERSION_MAJOR=4
ARG R_VERSION_MINOR=2
ARG R_VERSION_PATCH=2
ARG CONFIGURE_OPTIONS="--with-cairo --with-jpeglib --enable-R-shlib --with-blas --with-lapack"
RUN apt-get update && apt-get install -y -q --no-install-recommends \
            gfortran \
            git \
            git-annex \
            g++ \
            libreadline-dev \
            libx11-dev \
            libxt-dev \
            libpng-dev \
            libjpeg-dev \
            libcairo2-dev \   
            libssl-dev \ 
            libxml2-dev \
            libudunits2-dev \
            libgdal-dev \
            libbz2-dev \
            libzstd-dev \
            liblzma-dev \
            libpcre2-dev \
            locales \
            screen \
            texinfo \
            texlive \
            texlive-fonts-extra \
            vim \
            wget \
            xvfb \
            tcl8.6-dev \
            tk8.6-dev \
            cmake \
            curl \
            unzip \
            libcurl4-gnutls-dev \
            libgsl-dev \
            libcgal-dev \
            libglu1-mesa-dev \
            libglu1-mesa-dev \
            libtiff5-dev

RUN wget https://cran.rstudio.com/src/base/R-${R_VERSION_MAJOR}/R-${R_VERSION_MAJOR}.${R_VERSION_MINOR}.${R_VERSION_PATCH}.tar.gz \
            && tar zxvf R-${R_VERSION_MAJOR}.${R_VERSION_MINOR}.${R_VERSION_PATCH}.tar.gz \
            && rm R-${R_VERSION_MAJOR}.${R_VERSION_MINOR}.${R_VERSION_PATCH}.tar.gz \
            && cd /R-${R_VERSION_MAJOR}.${R_VERSION_MINOR}.${R_VERSION_PATCH} \
            && ./configure ${CONFIGURE_OPTIONS} \ 
            && make \
            && make install

RUN echo 'options(repos = c(CRAN = "https://cran.rstudio.com/"), download.file.method = "libcurl")' >> /usr/local/lib/R/etc/Rprofile.site \
            && Rscript -e "\
            install.packages(c('devtools', 'BiocManager', 'argparser', 'reticulate', 'rlist', 'oro.nifti', 'oro.dicom', 'fslr', 'WhiteStripe', 'matrixStats', 'R.matlab', 'abind', 'R.utils', 'RNifti', 'stapler', 'testthat', 'hexSticker', 'dplyr', 'oasis', 'fslr', 'plyr', 'misc3d', 'pixmap', 'colormap', 'ROCR', 'broom', 'broom.mixed', 'geepack', 'lme4', 'magrittr', 'neurobase', 'pbapply', 'purrr', 'readr', 'stringr', 'tibble', 'tidyr', 'voxel', 'forcats', 'gridExtra', 'RIA', 'Rfast', 'RJSONIO', 'pbmcapply', 'vesselr', 'tidyverse', 'caret', 'DMwR', 'openxlsx', 'acPCA', 'randomForest', 'RSpectra', 'doRNG', 'doParallel', 'doMC', 'doRedis', 'future.apply', 'badgecreatr'), repos='http://cran.rstudio.com'); \
            BiocManager::install(c('rhdf5', 'rhdf5filters', 'DelayedArray', 'DelayedMatrixStats', 'HDF5Array'));"

# Install latest version of cmake for ITKR (ANTsR dependency)
 RUN version=3.24 \
             && build=1 \
             && mkdir ~/temp \
             && cd ~/temp \
             && wget https://cmake.org/files/v$version/cmake-$version.$build.tar.gz \
             && tar -xzvf cmake-$version.$build.tar.gz \
             && cd cmake-$version.$build/ \
             && ./bootstrap \
             && make -j$(nproc) \
             && sudo make install \
             && alias cmake=$(which cmake) 
RUN git clone https://github.com/stnava/ITKR.git && R CMD INSTALL ITKR && rm -rf ITKR 
RUN git clone https://github.com/ANTsX/ANTsRCore.git && R CMD INSTALL ANTsRCore && rm -rf ANTsRCore
RUN git clone https://github.com/ANTsX/ANTsR.git && R CMD INSTALL ANTsR && rm -rf ANTsR

RUN Rscript -e "remotes::install_github('muschellij2/oasis', upgrade=TRUE)"
# Install packages from source
RUN git clone https://github.com/muschellij2/extrantsr.git && Rscript -e "install.packages('./extrantsr', repos = NULL, type='source')" && rm -rf extrantsr \
            && git clone https://github.com/muschellij2/NiftiArray.git && Rscript -e "install.packages('./NiftiArray', repos = NULL, type='source')" && rm -rf NiftiArray \
            && git clone https://github.com/avalcarcel9/mimosa.git && Rscript -e "install.packages('./mimosa', repos = NULL, type='source')" && rm -rf mimosa \
            && git clone https://github.com/muschellij2/malf.templates.git && Rscript -e "install.packages('./malf.templates', repos = NULL, type='source')" && rm -rf malf.templates 
# rtapas
RUN Rscript -e "install.packages(c('ggExtra', 'neuroim'))" \
            &&  git clone https://github.com/avalcarcel9/rtapas.git && Rscript -e "install.packages('./rtapas', repos = NULL, type='source')" && rm -rf rtapas
# lesiontools
RUN git clone https://github.com/jdwor/lesiontools.git && Rscript -e "install.packages('./lesiontools', repos = NULL, type='source')" && rm -rf lesiontools

# ANTsRNet
RUN Rscript -e "remotes::install_github('ANTsX/ANTsRNet', upgrade=TRUE)"

# kirby21 templates
RUN Rscript -e "install.packages(c('kirby21.base','kirby21.t1', 'kirby21.fmri'))" \
            -e "remotes::install_github('muschellij2/kirby21.flair')" 

# MNI template
RUN Rscript -e "remotes::install_github('Jfortin1/MNITemplate')"

# AnalyzeFMRI | Sublime
COPY AnalyzeFMRI_1.1-24.tar.gz .

RUN wget https://cran.r-project.org/src/contrib/fastICA_1.2-3.tar.gz && R CMD INSTALL fastICA_1.2-3.tar.gz && rm -rf fastICA_1.2-3.tar.gz \
    && wget https://cran.r-project.org/src/contrib/R.matlab_3.7.0.tar.gz && R CMD INSTALL R.matlab_3.7.0.tar.gz && rm -rf R.matlab_3.7.0.tar.gz \
    && R CMD INSTALL AnalyzeFMRI_1.1-24.tar.gz && rm -rf AnalyzeFMRI_1.1-24.tar.gz \
    && Rscript -e "remotes::install_github('emsweene/SuBLIME_package', upgrade=TRUE)"

# devtools
RUN apt-get update && apt-get install -y \
    libharfbuzz-dev \
    libfribidi-dev \
    && Rscript -e "install.packages('devtools')"

# hexSticker
RUN apt-get update && apt-get install -y \
    libmagick++-dev \
    && Rscript -e "install.packages('hexSticker')"

# colormap
RUN apt-get update && apt-get install -y \
    libv8-dev \
    && Rscript -e "install.packages('colormap')"

# DMwR
RUN Rscript -e "install.packages(c('xts', 'quantmod', 'zoo'))" \
    && wget https://cran.r-project.org/src/contrib/Archive/DMwR/DMwR_0.4.1.tar.gz && R CMD INSTALL DMwR_0.4.1.tar.gz && rm -rf DMwR_0.4.1.tar.gz

# doRedis
RUN apt-get update && apt-get install -y \
    libhiredis-dev \
    && Rscript -e "install.packages('doRedis')"

# badgecreatr
RUN Rscript -e "remotes::install_github('rmhogervorst/badgecreatr', upgrade=TRUE)"

WORKDIR /work
