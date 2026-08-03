FROM rockylinux:9 as builder

RUN dnf install -y 'dnf-command(config-manager)' && \
    dnf config-manager --enable crb && \
    dnf install -y git python3 mesa-libGL-devel ninja-build
WORKDIR /dials
COPY installer/bootstrap.py .
ENV PIP_ROOT_USER_ACTION=ignore
RUN python3 bootstrap.py --libtbx

# Copy to final image
FROM rockylinux:9
RUN dnf install -y glibc-locale-source
RUN localedef -i en_US -f UTF-8 en_US.UTF-8
RUN echo "LANG=\"en_US.UTF-8\"" > /etc/locale.conf
ENV LANG en_US.UTF-8

RUN mkdir /dials
COPY --from=builder /dials/conda_base /dials/conda_base
COPY --from=builder /dials/modules /dials/modules
COPY --from=builder /dials/build /dials/build
COPY --from=builder /dials/dials /dials
ENV PATH="/dials/conda_base/bin:/dials/build/bin:$PATH"
CMD ["dials.version"]
