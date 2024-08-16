FROM rockylinux:9 as builder

RUN dnf install -y 'dnf-command(config-manager)' && \
    dnf config-manager --enable crb && \
    dnf install -y git python3 mesa-libGL-devel ninja-build
WORKDIR /dials
COPY installer/bootstrap.py .
ENV PIP_ROOT_USER_ACTION=ignore
ENV CMAKE_GENERATOR=Ninja
RUN python3 bootstrap.py --cmake
RUN cmake --install build --config

# Copy to final image
FROM rockylinux:9
COPY ./docker-entrypoint.sh .
COPY --from=builder /dials/conda_base/dials/dials /dials/

RUN chmod 0755 /docker-entrypoint.sh

RUN dnf install -y glibc-locale-source 
RUN localedef -i en_US -f UTF-8 en_US.UTF-8
RUN echo "LANG=\"en_US.UTF-8\"" > /etc/locale.conf
ENV LANG en_US.UTF-8

ENTRYPOINT ["/docker-entrypoint.sh"]
CMD ["dials.version"]
