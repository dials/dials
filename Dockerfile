# Build dials
ARG BRANCH
FROM centos:7 as builder

RUN yum install -y git
WORKDIR /dials
COPY installer/bootstrap.py .
RUN python bootstrap.py \
    --branch "cctbx_project@${BRANCH}" \
    --branch "dials@${BRANCH}" \
    --branch "dxtbx@${BRANCH}" \
    --branch "cctbx_project@${BRANCH}"

# Copy to final image
FROM centos:7
COPY ./docker-entrypoint.sh .
COPY --from=builder /dials /dials

ENTRYPOINT ["/docker-entrypoint.sh"]
CMD ["dials.version"]