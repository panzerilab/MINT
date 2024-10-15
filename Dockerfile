# syntax=docker/dockerfile:1

FROM mathworks/matlab:r2023b
USER root
WORKDIR /app
COPY . .
RUN env | grep -i _PROXY
