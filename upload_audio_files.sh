extension="_spec-ratio_mle"
mle="MLE"

folder_local="Num mics/Num sources/N=2"
folder_remote="Num\ mics/Num\ sources/N=2"

# folder_local="Num sources"
# folder_remote="Num\ sources"

# Declare an array of string with type
#declare -a StringArray=("Ideal" "MCWF" "MLE" "MAP")
declare -a StringArray=("MLE")

 
# Iterate the string array using for loop
for val in ${StringArray[@]}; do
	if [ "$val" == "$mle" ]; then
		for sigma in 100
			do
				scp "./data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/${folder_local}/MLE/sigma=${sigma}/*${extension}.wav" orchi@ccrma-gate.stanford.edu:"~/Library/Web/Mic-bleed/Audio/${folder_remote}/MLE/sigma=${sigma}/"
		done
    else
    	scp "./data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/${folder_local}/${val}/*.wav" orchi@ccrma-gate.stanford.edu:"~/Library/Web/Mic-bleed/Audio/${folder_remote}/${val}/"
    fi
done