FROM ubuntu:18.04

MAINTAINER Luis Peralta

LABEL name="iwg"
LABEL description="Image to create an environment to run WEAVE IWG software"
LABEL version="0.01"

RUN apt-get update
RUN apt-get install -y python3-astropy python3-ipdb

COPY entrypoint.sh /

ENTRYPOINT ["/bin/bash", "/entrypoint.sh"]

CMD ["/bin/bash"]

