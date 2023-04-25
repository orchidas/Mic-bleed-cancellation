 #!/bin/bash

calib_type="spec-ratio"
method="map"
extension=${calib_type}_${method}
estimator="MAP" 
quartet_folder="26_TU_Berlin_Mozart_Quartett_Konzerthaus"


folder_remote="Mic\ Distance/Num\ Sources/N=2"
folder_local="Mic Distance/Num Sources/N=2"


# folder_remote="Num\ mics/Num\ sources/N=2"
# folder_local="Num mics/Num sources/N=2"

# folder_remote="RT60/Num\ sources/N=2"
# folder_local="RT60/Num sources/N=2"


# download saved results (but relevant files only)

scp orchi@cm-matlab.stanford.edu:"/scratch/orchi/src/data/${quartet_folder}/Audio/${folder_remote}/*${extension}*.mat" /
 "./data/${quartet_foder}/Audio/${folder_local}/"

for sigma in 0 1 100
	do
		scp orchi@cm-matlab.stanford.edu:"/scratch/orchi/src/data/${quartet_folder}/Audio/${folder_remote}/${estimator}/sigma=${sigma}/*${extension}.wav" /
		"./data/${quartet_folder}/Audio/${folder_local}/${estimator}/sigma=${sigma}/"
	done


# finished downloading
echo "Finished downloading results"
#quit remote server
