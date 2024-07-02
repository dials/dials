# Build dials
FROM rockylinux:8 as builder

RUN dnf install -y git python3.11 gcc gcc-c++ autoconf automake make gzip
WORKDIR /dials
COPY installer/bootstrap.py .
RUN python3 bootstrap.py

# Copy to final image
FROM rockylinux:8
COPY ./docker-entrypoint.sh .
COPY --from=builder /dials /dials
RUN chmod 0755 /docker-entrypoint.sh

RUN dnf install -y glibc-locale-source
RUN localedef -i en_US -f UTF-8 en_US.UTF-8
RUN echo "LANG=\"en_US.UTF-8\"" > /etc/locale.conf
ENV LANG en_US.UTF-8

ENTRYPOINT ["/docker-entrypoint.sh"]
CMD ["dials.version"]