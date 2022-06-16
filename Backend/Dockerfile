FROM node:16-alpine

WORKDIR /usr/src/app

COPY package.json package-lock.json ./

RUN npm install

RUN npm install -g typescript

COPY . .

RUN tsc

USER node

ENTRYPOINT node -r source-map-support/register --trace-uncaught build/server.js