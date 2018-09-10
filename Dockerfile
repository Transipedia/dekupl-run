FROM python:3.6-stretch

LABEL maintainer="Dimitri Larue <dimitri.larue@seq.one>"

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    less \
    locales \
    wget \
    ca-certificates \
    fonts-texgyre \
    jellyfish \
    pigz \
    cmake \
    perl \
  && rm -rf /var/lib/apt/lists/*

# R
ENV R_BASE_VERSION 3.3.3


## Now install R and littler, and create a link for littler in /usr/local/bin
## Also set a default CRAN repo, and make sure littler knows about it too
## Also install stringr to make dococt install (from source) easier
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    littler \
    r-cran-littler \
    r-cran-stringr \
    r-base=${R_BASE_VERSION}-* \
    r-base-dev=${R_BASE_VERSION}-* \
    r-recommended=${R_BASE_VERSION}-* \
      && echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"))' >> /etc/R/Rprofile.site \
      && echo 'source("/etc/R/Rprofile.site")' >> /etc/littler.r \
  && ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
  && ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
  && ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
  && ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
  && install.r docopt \
  && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
&& rm -rf /var/lib/apt/lists/*

WORKDIR /dekupl
COPY install_r_packages.R .
RUN Rscript install_r_packages.R

# Python
RUN pip install "snakemake<5.0.0"

COPY bin bin
COPY share share

RUN cd share/computeNF \
  && rm -f /dekupl/bin/computeNF \
  && make \
  && cp computeNF /dekupl/bin/ \
  && cd ../../share/joinCounts \
  && rm -f /dekupl/bin/joinCounts \
  && make \
  && cp joinCounts /dekupl/bin/ \
  && cd ../../share/mergeTags \
  && rm -f /dekupl/bin/mergeTags \
  && make \
  && cp mergeTags /dekupl/bin/ \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
  libboost-all-dev \
  && cd ../../share/TtestFilter \
  && rm -f /dekupl/bin/TtestFilter \
  && make \
  && cp TtestFilter /dekupl/bin/ \
  && rm -rf /var/lib/apt/lists/* \
  && cd ../../share \
  && rm -f /dekupl/bin/kallisto \
  && wget https://github.com/pachterlab/kallisto/releases/download/v0.43.0/kallisto_linux-v0.43.0.tar.gz -O kallisto.tar.gz \
  && tar -xzf kallisto.tar.gz \
  && cp kallisto_linux-v0.43.0/kallisto /dekupl/bin
RUN rm -rf share

COPY Snakefile .

ENTRYPOINT [ "snakemake", "-s", "/dekupl/Snakefile" ]