FROM dasteck/cglass:separation
WORKDIR /build
RUN git clone --recursive --single-branch --branch separation_paper https://github.com/Betterton-Lab/C-GLASS.git . &&\
    ./install.sh -otI
RUN rm -rf /build
WORKDIR /mnt
ENTRYPOINT []
