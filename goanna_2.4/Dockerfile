## Dockerfile
FROM ubuntu:22.04
MAINTAINER Amanda Cooksey
LABEL Description="AgBase GOanna"

ENV DEBIAN_FRONTEND=noninteractive

# Install all the updates and download dependencies
RUN apt-get update && \
    apt-get install -y \
    git \
    wget \
    python3 \
    python3-scipy \
    python3-sklearn \
    python3-numpy \
    python3-biopython \
    locales \
    bzip2

RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && locale-gen

ENV LANG="en_US.UTF-8" LANGUAGE="en_US:en" LC_ALL="en_US.UTF-8"

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_25.7.0-2-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh


# give write permissions to conda folder
RUN chmod 777 -R /opt/conda/

ENV PATH=$PATH:/opt/conda/bin

RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

RUN conda config --add channels bioconda

RUN conda upgrade conda

# add blast
RUN conda install -c conda-forge -c bioconda blast taxonkit && conda install pandas
#RUN conda install -c conda-forge -c bioconda taxonkit
#RUN conda install pandas

ENV PATH /masterscript.sh/:$PATH

# add masterscript--i added this
ADD masterscript.sh splitB.pl cyverse_blast2GO.pl taxon_constrain.py /usr/bin/

#ADD splitB.pl /usr/bin

#ADD cyverse_blast2GO.pl /usr/bin

#ADD taxon_constrain.py /usr/bin

# Change the permissions and the path for the wrapper script
RUN chmod +x /usr/bin/masterscript.sh /usr/bin/taxon_constrain.py

RUN mkdir -p /work-dir /agbase_database /go_info

RUN chmod -R a+rwx /work-dir /agbase_database /go_info

ADD GO0032991_and_children.json /usr

# Entrypoint
ENTRYPOINT ["/usr/bin/masterscript.sh"]

# Add path to working directory
WORKDIR /work-dir
