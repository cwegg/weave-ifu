FROM ubuntu:18.04

MAINTAINER Luis Peralta

LABEL name="iwg"
LABEL description="Image to create an environment to run WEAVE IWG software"
LABEL version="0.03"

RUN apt-get update
RUN apt-get install -y python3-ipdb python3-astropy astropy-utils

COPY entrypoint.sh /

ENTRYPOINT ["/bin/bash", "/entrypoint.sh"]

CMD ["/bin/bash"]

