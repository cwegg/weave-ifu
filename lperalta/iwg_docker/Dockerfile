FROM ubuntu:20.04

MAINTAINER Luis Peralta

LABEL name="iwg"
LABEL description="Image to create an environment to run WEAVE IWG software"
LABEL version="0.04"

RUN apt-get update
RUN apt-get install -y python3-ipdb python3-astropy astropy-utils python3-astropy-healpix python3-matplotlib

COPY entrypoint.sh /

ENTRYPOINT ["/bin/bash", "/entrypoint.sh"]

CMD ["/bin/bash"]

