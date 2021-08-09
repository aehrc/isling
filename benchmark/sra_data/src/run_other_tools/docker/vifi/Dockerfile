#  $ docker build . -t szsctt/vifi:latest -t szsctt/vifi:1
#  $ docker run --rm -it szsctt/vifi:latest /bin/bash
#  $ docker push szsctt/vifi:latest
#  $ docker push szsctt/vifi:1

# cd ~/Documents/Projects/viInt/experiments/expt6_simulations/vifi
# docker run -v$(pwd):/home --rm szsctt/vifi:1

FROM namphuon/vifi

COPY ViFi/scripts /home/ViFi/scripts

ENV VIFI_DIR=/home/ViFi
