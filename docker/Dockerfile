FROM jeffmm/simcore_base:latest
WORKDIR /build
RUN git clone --recursive --single-branch --branch master https://github.com/Betterton-Lab/C-GLASS.git . &&\
    ./install.sh -otI
RUN rm -rf /build
WORKDIR /mnt
ENTRYPOINT []
