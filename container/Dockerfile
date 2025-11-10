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

RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip -O /tmp/fastqc.zip && \ 
                unzip /tmp/fastqc.zip -d /opt/ && \
                rm /tmp/fastqc.zip && \
                chmod 777 /opt/FastQC/fastqc

ENV PATH="/opt/FastQC/:${PATH}"

RUN Rscript /app/install_r_packages && \
    pip install -r /app/requirements.txt
