FROM amanacevice/pandas

RUN apt-get install -y net-tools
RUN apt-get update -qq && apt-get -y install wget


RUN pip3 install pcst_fast numpy matplotlib python-igraph goenrich leidenalg

RUN python setup.py install

##now get PPI network and create
RUN mkdir data
RUN wget -O data/ppi_net.txt.gz https://stringdb-static.org/download/protein.links.detailed.v11.0/9606.protein.links.detailed.v11.0.txt.gz
RUN wget -O data/ppi_map.txt.gz https://stringdb-static.org/download/protein.info.v11.0/9606.protein.info.v11.0.txt.gz
RUN python ./bin/graphTools.py --graphFile data/ppi_net.txt.gz --graphSource=string --nodeMapping=ppi_map.txt.gz --dest=./data

VOLUME ['/tmp']
ENTRYPOINT ["python", "/bin/createHyphaFromProts.py"]
