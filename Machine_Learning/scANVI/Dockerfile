FROM python:3.9-slim

ENV APP_HOME /app
WORKDIR $APP_HOME
ADD main.py .
ADD utils.py .
ADD requirements.txt .

# Install production dependencies.
RUN pip install --no-cache-dir -r requirements.txt