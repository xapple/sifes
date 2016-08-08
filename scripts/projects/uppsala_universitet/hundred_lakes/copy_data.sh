# From milou to warwick #
rsync -av --progress /proj/b2011035/INBOX/160321_M00485_0262_000000000-AMK6/ warwick-node:/home/lucas/SIFES/raw/projects/uppsala_universitet/hundred_lakes/

# Remove write permission #
find ~/SIFES/raw/projects/envonautics/thaw_ponds -type d -print0 | xargs -0 chmod u=rx,g=,o=
find ~/SIFES/raw/projects/envonautics/thaw_ponds -type f -print0 | xargs -0 chmod u=r,g=,o=
