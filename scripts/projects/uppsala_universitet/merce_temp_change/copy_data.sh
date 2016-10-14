# From milou to warwick (proj b2014083) #
rsync -av --progress --exclude "Plots" /proj/b2014083/INBOX/150612_M00485_0205_000000000-AFJ7P/ warwick-node:/home/lucas/SIFES/raw/projects/uppsala_universitet/merce_temp_change/

# Remove write permission #
find /home/lucas/SIFES/raw/projects/uppsala_universitet/merce_temp_change/ -type d -print0 | xargs -0 chmod u=rx,g=,o=
find /home/lucas/SIFES/raw/projects/uppsala_universitet/merce_temp_change/ -type f -print0 | xargs -0 chmod u=r,g=,o=

# Generate JSON files (was originally run016) #
~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/uppsala_universitet/merce_temp_change/metadata.xlsx

# Copy back from storage 1 to storage 2 #
rsync -av --progress /home/lucas/storage1/sifes/raw/projects/uppsala_universitet/merce_temp_change/ /home/lucas/storage2/sifes/raw/projects/uppsala_universitet/merce_temp_change/
rsync -av --progress /home/lucas/storage1/sifes/views/projects/uppsala_universitet/merce_temp_change/ /home/lucas/storage2/sifes/views/projects/uppsala_universitet/merce_temp_change/
rsync -av --progress /home/lucas/storage1/sifes/views/samples/uppsala_universitet/merce_temp_change/ /home/lucas/storage2/sifes/views/samples/uppsala_universitet/merce_temp_change/