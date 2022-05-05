FROM minio/mc:latest

ENV MINIO_HOST=http://minio:9000
ENV MINIO_USER=minioadmin
ENV MINIO_PASSWORD=minioadmin


ENTRYPOINT \
mc config host add devminio ${MINIO_HOST} ${MINIO_USER} ${MINIO_PASSWORD}; \
mc mb devminio/minio-bucket; \
mc mb devminio/minio-picture-bucket; \
mc policy set download devminio/minio-picture-bucket; \
exit 0;
