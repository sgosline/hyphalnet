FROM python:3.7

#RUN apt-get install -y net-tools\\\\\\\\\\\\\\
#RUN apk update -qq && \
#        apk add wget && \
#        apk add libtool && \
#        apk add m4 && \
#        apk add autoconf &&\
#        apk add automake && \
#        apk add patch &&\
#        apk add gcc &&\
#        apk add g++
RUN apt-get update -qq
##        && apt-get install -y wget

COPY . hyphalnet
WORKDIR hyphalnet
RUN ls -la *
RUN python setup.py install

##now get PPI network and create
RUN mkdir data
RUN wget -O data/ppi_net.txt.gz https://stringdb-static.org/download/protein.links.detailed.v11.0/9606.protein.links.detailed.v11.0.txt.gz
RUN wget -O data/ppi_map.txt.gz https://stringdb-static.org/download/protein.info.v11.0/9606.protein.info.v11.0.txt.gz
RUN python bin/graphTools.py --graphFile data/ppi_net.txt.gz --graphSource=string --nodeMapping data/ppi_map.txt.gz --dest=data

#RUN mv /data/* data && \
#        rmdir /data/

VOLUME ['/tmp']
ENTRYPOINT ["python", "bin/createHyphaFromProts.py"]
