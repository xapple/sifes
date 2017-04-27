#!/usr/bin/env bash

# To be run from enda #
#for m in proj.cluster.locations_maps: print m.path
#for s in proj: print s.graphs.location_map.path

# To be run from warwick #
chmod +x ~/deploy/sifes/scripts/projects/unige/desalt/copy_maps.sh

mkdir -p ~/SIFES/views/projects/unige/desalt/cluster/desalt/graphs/

mkdir -p ~/SIFES/views/samples/unige/desalt/as1a/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/as1c/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/as2a/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/as2b/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/as2c/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/as5a/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/as5b/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/as5c/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/vmso1a/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/vmso1b/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/vmso1c/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/vmso2a/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/vmso2b/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/vmso2c/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/vmso3a/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/vmso3b/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/vmso3c/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/had1b/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/had1c/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/had2a/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/had2b/graphs
mkdir -p ~/SIFES/views/samples/unige/desalt/had3a/graphs

rsync -avz edna:/home/sinclair/SIFES/views/projects/unige/desalt/cluster/desalt/graphs/location_map_ashqelon.png ~/SIFES/views/projects/unige/desalt/cluster/desalt/graphs/location_map_ashqelon.png
rsync -avz edna:/home/sinclair/SIFES/views/projects/unige/desalt/cluster/desalt/graphs/location_map_soreq.png    ~/SIFES/views/projects/unige/desalt/cluster/desalt/graphs/location_map_soreq.png
rsync -avz edna:/home/sinclair/SIFES/views/projects/unige/desalt/cluster/desalt/graphs/location_map_hadera.png   ~/SIFES/views/projects/unige/desalt/cluster/desalt/graphs/location_map_hadera.png

rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/as1a/graphs/location_map.png   ~/SIFES/views/samples/unige/desalt/as1a/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/as1c/graphs/location_map.png   ~/SIFES/views/samples/unige/desalt/as1c/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/as2a/graphs/location_map.png   ~/SIFES/views/samples/unige/desalt/as2a/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/as2b/graphs/location_map.png   ~/SIFES/views/samples/unige/desalt/as2b/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/as2c/graphs/location_map.png   ~/SIFES/views/samples/unige/desalt/as2c/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/as5a/graphs/location_map.png   ~/SIFES/views/samples/unige/desalt/as5a/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/as5b/graphs/location_map.png   ~/SIFES/views/samples/unige/desalt/as5b/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/as5c/graphs/location_map.png   ~/SIFES/views/samples/unige/desalt/as5c/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/vmso1a/graphs/location_map.png ~/SIFES/views/samples/unige/desalt/vmso1a/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/vmso1b/graphs/location_map.png ~/SIFES/views/samples/unige/desalt/vmso1b/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/vmso1c/graphs/location_map.png ~/SIFES/views/samples/unige/desalt/vmso1c/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/vmso2a/graphs/location_map.png ~/SIFES/views/samples/unige/desalt/vmso2a/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/vmso2b/graphs/location_map.png ~/SIFES/views/samples/unige/desalt/vmso2b/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/vmso2c/graphs/location_map.png ~/SIFES/views/samples/unige/desalt/vmso2c/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/vmso3a/graphs/location_map.png ~/SIFES/views/samples/unige/desalt/vmso3a/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/vmso3b/graphs/location_map.png ~/SIFES/views/samples/unige/desalt/vmso3b/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/vmso3c/graphs/location_map.png ~/SIFES/views/samples/unige/desalt/vmso3c/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/had1b/graphs/location_map.png  ~/SIFES/views/samples/unige/desalt/had1b/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/had1c/graphs/location_map.png  ~/SIFES/views/samples/unige/desalt/had1c/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/had2a/graphs/location_map.png  ~/SIFES/views/samples/unige/desalt/had2a/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/had2b/graphs/location_map.png  ~/SIFES/views/samples/unige/desalt/had2b/graphs/location_map.png
rsync -avz edna:/home/sinclair/SIFES/views/samples/unige/desalt/had3a/graphs/location_map.png  ~/SIFES/views/samples/unige/desalt/had3a/graphs/location_map.png

