#!/bin/sh

extension="_spec-ratio_map"
estimator="MAP"


# download saved results (but relevant files only)

for sigma in 0 1 100
	do
		scp orchi@cm-matlab.stanford.edu:"/scratch/orchi/src/data/Drums_Noah/Saved\ files/separated_drums${extension}_sigma=${sigma}.mat" "./data/Drums_Noah/Saved files/"
		scp orchi@cm-matlab.stanford.edu:"/scratch/orchi/src/data/Drums_Noah/Saved\ files/${estimator}/sigma=${sigma}/*${extension}.wav" "./data/Drums_Noah/Saved files/${estimator}/sigma=${sigma}/"
	done

# finished downloading
echo "Finished downloading results"
#quit remote server
