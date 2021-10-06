 #!/bin/bash

calib_type="spec-ratio"
method="map"
extension=${calib_type}_${method}
estimator="MAP" 


# folder_remote="Mic\ Distance/Num\ Sources/N=2"
# folder_local="Mic Distance/Num Sources/N=2"


# folder_remote="Num\ mics/Num\ sources/N=2"
# folder_local="Num mics/Num sources/N=2"

folder_remote="RT60/Num\ sources/N=2"
folder_local="RT60/Num sources/N=2"


# download saved results (but relevant files only)
# scp -r orchi@cm-matlab.stanford.edu:"/scratch/orchi/src/data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/Mic\ Distance/Num\ Sources/N=2/*spec-ratio_mle.mat" ./data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/Mic\ Distance/Num\ Sources/
scp orchi@cm-matlab.stanford.edu:"/scratch/orchi/src/data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/${folder_remote}/*${extension}*.mat" "./data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/${folder_local}/"

for sigma in 0 1 100
	do
		scp orchi@cm-matlab.stanford.edu:"/scratch/orchi/src/data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/${folder_remote}/${estimator}/sigma=${sigma}/*${extension}.wav" "./data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/${folder_local}/${estimator}/sigma=${sigma}/"
	done


# finished downloading
echo "Finished downloading results"
#quit remote server
