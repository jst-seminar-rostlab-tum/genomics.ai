FROM minio/minio:latest

RUN mkdir /data
RUN mkdir /data/minio-bucket
RUN mkdir /data/minio-picture-bucket
COPY default-user.png /data/minio-picture-bucket/user.png

CMD ["server", "--console-address", ":9001", "/data"]
