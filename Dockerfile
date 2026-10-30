FROM python:3.12.16-slim-bookworm@sha256:623674a14ef93d052bb9b0de4dbf5f28eb578a35263838941cc8aa43213c1487

ENV PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /processor

# Pinned Python dependencies (wheels only), plus the Lambda runtime interface client for running
# as a Lambda function. No scanpy: the counts are read with h5py and anndata's element reader.
COPY requirements.txt /processor/requirements.txt
RUN pip install --only-binary=:all: -r /processor/requirements.txt awslambdaric==2.2.1 simplejson==4.2.0 \
    && find /usr/local/lib/python3.12/site-packages -name '*.pyc' -delete \
    && find /usr/local/lib/python3.12/site-packages -type d -name tests -prune -exec rm -rf {} +

COPY processor/ /processor/processor/
COPY VERSION /processor/VERSION

ENV PYTHONPATH="/processor" \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib

# Build matplotlib's font cache now (as the run-time user would), so a run needs no writes outside /tmp.
RUN python -c "import matplotlib.pyplot, pydeseq2.dds, anndata" && chmod -R a+rwX /tmp/matplotlib

# Runs as uid 1000 on the platform; nothing here needs root.
USER 1000

# Runs main() once as a container, or serves processor.handler as a Lambda function (see processor/entry.py).
ENTRYPOINT ["python3", "-m", "processor.entry"]
