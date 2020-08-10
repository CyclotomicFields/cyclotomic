FROM debian:latest
COPY scripts /cyclotomic_scripts
RUN apt-get update && \
    apt-get install -y sudo wget build-essential automake autoconf curl
RUN cd /cyclotomic_scripts && ./install_flint.sh
RUN cd /cyclotomic_scripts && ./install_antic.sh
RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
RUN rustup default nightly
