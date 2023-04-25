#!/bin/sh

extension="_spec-ratio_map"
estimator="MAP"
drums_folder="Drums_Noah"


# download saved results (but relevant files only)

for sigma in 0 1 100
	do
		scp orchi@cm-matlab.stanford.edu:"/scratch/orchi/src/data/${drum_folder}/Saved\ files/separated_drums${extension}_sigma=${sigma}.mat" \
		"./data/${drum_folder}/Saved files/"

		scp orchi@cm-matlab.stanford.edu:"/scratch/orchi/src/data/${drum_folder}/Saved\ files/${estimator}/sigma=${sigma}/*${extension}.wav" \
		 		"./data/${drum_folder}/Saved files/${estimator}/sigma=${sigma}/"
	done

# finished downloading
echo "Finished downloading results"
#quit remote server
