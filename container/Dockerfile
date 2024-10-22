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
				texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texlive-latex-extra lmodern

RUN Rscript /app/install_r_packages && \
    pip install -r /app/requirements.txt
