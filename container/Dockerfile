FROM python:3.11-bookworm

RUN mkdir /app

COPY app/* /app/

ENV R_BASE_VERSION=4.2.2

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
                r-base=4.2.2.20221110-2 \
                r-base-dev=4.2.2.20221110-2 \
                r-base-core=4.2.2.20221110-2 \
				vim \
				r-cran-devtools \
                default-jre \
				texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texlive-latex-extra lmodern

RUN apt-get -qq update && apt-get -y upgrade && \
                apt install -y wget libfindbin-libs-perl software-properties-common unzip

RUN apt-get update && apt-get install -y \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN curl -sSL "https://github.com/peak/s5cmd/releases/download/v2.3.0/s5cmd_2.3.0_Linux-64bit.tar.gz" | tar -xz -C /usr/local/bin s5cmd

RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip -O /tmp/fastqc.zip && \ 
                unzip /tmp/fastqc.zip -d /opt/ && \
                rm /tmp/fastqc.zip && \
                chmod 777 /opt/FastQC/fastqc

RUN apt-get install -y curl unzip \
    && curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
    && unzip awscliv2.zip \
    && ./aws/install \
    && rm -rf awscliv2.zip aws

ENV PATH="/opt/FastQC/:${PATH}"

RUN Rscript /app/install_r_packages && \
    pip install -r /app/requirements.txt
