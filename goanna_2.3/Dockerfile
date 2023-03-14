## Dockerfile
FROM ubuntu:22.04
LABEL Description="AgBase GOanna version 2.3"

# Install all the updates and download dependencies
RUN apt-get update && \
    apt-get install -y \
    locales \
    git \
    wget \
    bzip2

#RUN locale-gen en_US.UTF-8
# Set the locale
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && \
    locale-gen
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_23.1.0-1-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh


# give write permissions to conda folder
RUN chmod 777 -R /opt/conda/

ENV PATH=$PATH:/opt/conda/bin

RUN conda config --add channels bioconda

RUN conda upgrade conda

# ADD BLAST 2.13.0
RUN conda install -c conda-forge -c bioconda blast==2.13.0

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
