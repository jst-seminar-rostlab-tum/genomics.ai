FROM minio/minio:latest

RUN mkdir /data
RUN mkdir /data/minio-bucket

CMD ["server", "--console-address", ":9001", "/data"]