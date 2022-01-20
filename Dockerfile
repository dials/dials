# Build dials
FROM centos:7 as builder

RUN yum install -y git
WORKDIR /dials
COPY installer/bootstrap.py .
RUN python bootstrap.py

# Copy to final image
FROM centos:7
COPY ./docker-entrypoint.sh .
COPY --from=builder /dials /dials
RUN chmod 0755 /docker-entrypoint.sh

ENTRYPOINT ["/docker-entrypoint.sh"]
CMD ["dials.version"]