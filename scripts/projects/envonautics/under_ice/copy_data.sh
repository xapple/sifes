# From uppmax #
rsync -av --progress taito4:/wrk/eiler/BaseSpace/Experiment_1/160215-EHa_Stromquist_AmplSeq1run-29109084/ ~/SIFES/raw/projects/micans/micans_v6/experiment_1/

rsync -av --progress taito4:/wrk/eiler/BaseSpace/Experiment_2/160215-EHa_Stromquist_AmplSeq1run-29109084/ ~/SIFES/raw/projects/micans/micans_v6/experiment_2/

# Remove write permission #
find ~/SIFES/raw/projects/micans/micans_v6 -type d -print0 | xargs -0 chmod u=rx,g=,o=
find ~/SIFES/raw/projects/micans/micans_v6 -type f -print0 | xargs -0 chmod u=r,g=,o=
