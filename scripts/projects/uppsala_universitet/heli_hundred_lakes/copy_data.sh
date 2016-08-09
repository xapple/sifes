# From milou to warwick (proj b2011035 orig) #
rsync -av --progress  --exclude "Plots" /proj/b2011035/INBOX/160321_M00485_0262_000000000-AMK6/ warwick-node:/home/lucas/SIFES/raw/projects/uppsala_universitet/hundred_lakes/

# From milou to warwick (proj copied) #
rsync -av --progress --exclude "Plots" /proj/b2015084/INBOX/160321_M00485_0262_000000000-AMK63/ warwick-node:/home/lucas/SIFES/raw/projects/uppsala_universitet/hundred_lakes/

# Remove write permission #
find /home/lucas/SIFES/raw/projects/uppsala_universitet/hundred_lakes/ -type d -print0 | xargs -0 chmod u=rx,g=,o=
find /home/lucas/SIFES/raw/projects/uppsala_universitet/hundred_lakes/ -type f -print0 | xargs -0 chmod u=r,g=,o=

# Generate JSON files #
~/repos/sifes/metadata/excel_to_json.py ~/repos/sifes/metadata/excel/projects/uppsala_universitet/heli_hundred_lakes/metadata.xlsx