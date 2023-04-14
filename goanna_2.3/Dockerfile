## Dockerfile
FROM ubuntu:16.04
MAINTAINER Amanda Cooksey
LABEL Description="AgBase GOanna"

# Install all the updates and download dependencies
RUN apt-get update && \
    apt-get install -y \
    git \
    wget \
    bzip2

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.0.5-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh


# give write permissions to conda folder
RUN chmod 777 -R /opt/conda/

ENV PATH=$PATH:/opt/conda/bin

RUN conda config --add channels bioconda

RUN conda upgrade conda

<<<<<<< HEAD
# add blast 2.7.1
RUN conda install -c conda-forge -c bioconda blast==2.7.1
=======
# ADD BLAST 2.13.0
RUN conda install -c conda-forge -c bioconda blast==2.13.0
>>>>>>> origin/master

ENV PATH /masterscript.sh/:$PATH

# add masterscript--i added this
ADD masterscript.sh /usr/bin

ADD splitB.pl /usr/bin

ADD cyverse_blast2GO.pl /usr/bin

# Change the permissions and the path for the wrapper script
RUN chmod +x /usr/bin/masterscript.sh

RUN mkdir /work-dir
RUN mkdir /agbase_database
RUN mkdir /go_info


ADD GO0032991_and_children.json /usr

# Entrypoint
ENTRYPOINT ["/usr/bin/masterscript.sh"]

# Add path to working directory
WORKDIR /work-dir
