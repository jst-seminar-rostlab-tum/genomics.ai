FROM python:3.9-slim

RUN apt-get -y update && \
    apt-get -y install gcc&& \
    apt-get clean

ENV APP_HOME /app
WORKDIR $APP_HOME

ADD requirements.txt .

RUN pip install -r requirements.txt --no-cache

COPY . .
ENV WORKERS 1
ENV THREADS 1
ENV PORT 8080

CMD exec gunicorn --bind :$PORT --workers $WORKERS --threads $THREADS --timeout 0 api:app