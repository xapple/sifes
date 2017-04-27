#!/usr/bin/env bash

# To be run from enda #
#for m in proj.cluster.locations_maps: print m.path
#for s in proj: print s.graphs.location_map.path

# To be run from warwick #
chmod +x ~/deploy/sifes/scripts/projects/unige/foram/copy_maps.sh

mkdir -p ~/SIFES/views/projects/unige/foram/cluster/foram/graphs/

mkdir -p ~/SIFES/views/samples/unige/foram/as1a/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/as1b/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/as1c/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/as2a/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/as2b/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/as2c/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/as5b/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/as5c/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/vmso1a/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/vmso1b/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/vmso1c/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/vmso2a/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/vmso2c/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/vmso3a/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/vmso3c/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/had1b/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/had2a/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/had2b/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/had2c/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/had3a/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/had3b/graphs
mkdir -p ~/SIFES/views/samples/unige/foram/had3c/graphs

rsync -avz edna:/home/sinclair/SIFES/views/projects/unige/foram/cluster/foram/graphs/location_map_ashqelon.png ~/SIFES/views/projects/unige/foram/cluster/foram/graphs/location_map_ashqelon.png
rsync -avz edna:/home/sinclair/SIFES/views/projects/unige/foram/cluster/foram/graphs/location_map_soreq.png    ~/SIFES/views/projects/unige/foram/cluster/foram/graphs/location_map_soreq.png
rsync -avz edna:/home/sinclair/SIFES/views/projects/unige/foram/cluster/foram/graphs/location_map_hadera.png   ~/SIFES/views/projects/unige/foram/cluster/foram/graphs/location_map_hadera.png

rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/as1a/graphs/location_map.png   ~/SIFES/views/samples/unige/foram/as1a/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/as1b/graphs/location_map.png   ~/SIFES/views/samples/unige/foram/as1b/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/as1c/graphs/location_map.png   ~/SIFES/views/samples/unige/foram/as1c/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/as2a/graphs/location_map.png   ~/SIFES/views/samples/unige/foram/as2a/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/as2b/graphs/location_map.png   ~/SIFES/views/samples/unige/foram/as2b/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/as2c/graphs/location_map.png   ~/SIFES/views/samples/unige/foram/as2c/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/as5b/graphs/location_map.png   ~/SIFES/views/samples/unige/foram/as5b/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/as5c/graphs/location_map.png   ~/SIFES/views/samples/unige/foram/as5c/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/vmso1a/graphs/location_map.png ~/SIFES/views/samples/unige/foram/vmso1a/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/vmso1b/graphs/location_map.png ~/SIFES/views/samples/unige/foram/vmso1b/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/vmso1c/graphs/location_map.png ~/SIFES/views/samples/unige/foram/vmso1c/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/vmso2a/graphs/location_map.png ~/SIFES/views/samples/unige/foram/vmso2a/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/vmso2c/graphs/location_map.png ~/SIFES/views/samples/unige/foram/vmso2c/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/vmso3a/graphs/location_map.png ~/SIFES/views/samples/unige/foram/vmso3a/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/vmso3c/graphs/location_map.png ~/SIFES/views/samples/unige/foram/vmso3c/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/had1b/graphs/location_map.png  ~/SIFES/views/samples/unige/foram/had1b/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/had2a/graphs/location_map.png  ~/SIFES/views/samples/unige/foram/had2a/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/had2b/graphs/location_map.png  ~/SIFES/views/samples/unige/foram/had2b/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/had2c/graphs/location_map.png  ~/SIFES/views/samples/unige/foram/had2c/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/had3a/graphs/location_map.png  ~/SIFES/views/samples/unige/foram/had3a/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/had3b/graphs/location_map.png  ~/SIFES/views/samples/unige/foram/had3b/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/foram/had3c/graphs/location_map.png  ~/SIFES/views/samples/unige/foram/had3c/graphs/location_map.png

