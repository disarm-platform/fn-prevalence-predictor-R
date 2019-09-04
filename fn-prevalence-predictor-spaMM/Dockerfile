FROM rocker/geospatial

RUN wget https://github.com/openfaas/faas/releases/download/0.9.14/fwatchdog -O /usr/bin/fwatchdog
RUN chmod +x /usr/bin/fwatchdog

RUN mkdir -p /home/app/function

WORKDIR /app

COPY function/install_packages.R function/
RUN Rscript function/install_packages.R

COPY main.R .
COPY function function

# Populate example here - i.e. "cat", "sha512sum" or "node index.js"
ENV fprocess="Rscript main.R"
# Set to true to see request in function logs
# ENV write_debug="false"

EXPOSE 8080

HEALTHCHECK --interval=3s CMD [ -e /tmp/.lock ] || exit 1
CMD [ "fwatchdog" ]
