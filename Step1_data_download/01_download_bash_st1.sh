#!/bin/bash

#code to download the full extent of our gbif dataset
#cd /scratch/kricks/GBIF_store/processed_Rhizobia/
#note this can be downloaded through nohup which is way faster
wget https://api.gbif.org/v1/occurrence/download/request/0011505-260806074905277.zip.zip
#nohup wget https://api.gbif.org/v1/occurrence/download/request/0017640-260721160103020.zip &

#unzip and extract enormous file
mkdir -p gbif_extracted
unzip -o 0011505-260806074905277.zip -d gbif_extracted
ls -lh gbif_extracted