FROM registry.gitlab.com/transipedia/dekupl-run:base

WORKDIR /dekupl
COPY bin bin
RUN ln -s /usr/local/bin/computeNF /dekupl/bin/computeNF \
  && ln -s /usr/local/bin/joinCounts /dekupl/bin/joinCounts \
  && ln -s /usr/local/bin/mergeTags /dekupl/bin/mergeTags \
  && ln -s /usr/local/bin/TtestFilter /dekupl/bin/TtestFilter \
  && ln -s /usr/local/bin/kallisto /dekupl/bin/kallisto

COPY config.json .
COPY Snakefile .

ENTRYPOINT [ "snakemake", "-s", "/dekupl/Snakefile" ]
