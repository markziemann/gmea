FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Update apt-get
RUN apt-get update \
        && apt-get upgrade -y \
        && apt-get install -y nano git  libncurses-dev \
        ## Install the python package magic wormhole to send files
        && pip install magic-wormhole           \
        ## Remove packages in '/var/cache/' and 'var/lib'
        ## to remove side-effects of apt-get update
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/*

# Install CRAN packages
RUN Rscript -e 'install.packages(c("beeswarm","devtools","dplyr","eulerr","forestplot","gplots","gridExtra","HGNChelper","kableExtra","parallel","plyr","png","psych","qqman","RCircos","RColorBrewer","reshape2","R.utils","stringi","tictoc"))'

# Install bioconductor packages
RUN Rscript -e 'BiocManager::install(c("DESeq2","DMRcate","DMRcatedata","ENmix","GEOquery","IlluminaHumanMethylation450kanno.ilmn12.hg19","IlluminaHumanMethylation450kmanifest","IlluminaHumanMethylationEPICanno.ilm10b4.hg19","limma","minfi","missMethyl","org.Hs.eg.db","topconfects"))'

# install mitch separately so that the correct aggregation function is used
RUN Rscript -e 'devtools::install_github("markziemann/mitch")'

# get a clone of the codes
RUN git clone https://github.com/markziemann/gmea.git

# Set the container working directory
ENV DIRPATH /gmea
WORKDIR $DIRPATH

