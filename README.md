<h2>Microphone cross-talk cancellation in ensemble recordings</h2>

<h3>Manual</h3>

<p> This repository contains code related to my PhD thesis at CCRMA. See <a href="my code/debleed.m">debleed.m</a> for how to use it. For sound examples, visit <a href = "https://ccrma.stanford.edu/~orchi/Mic-bleed/sound_examples_quartet.html">string quartet</a> and <a href = "https://ccrma.stanford.edu/~orchi/Mic-bleed/sound_examples_drums.html">drums</a>.</p>

<h3>Abstract</h3>

<p>While recording an ensemble of musicians, microphone cross-talk, or “bleed”, is considered a nuisance by audio engineers. When two microphones pick up the same signal with a time delay, comb filtering artifacts are present. In the close-miking technique, microphones are placed at a distance of 5 - 50 cm from the sound source. While this attenuates intereferences to some extent, cross-talk can become significant in such cases due to the effect of
nearby strong reflective surfaces. The complexity lies in the fact that there is usually an arbitrarynumber and distribution of instruments and microphones, and results are influenced by the room acoustics of the studio where the ensemble is recorded.</p>

<p>In this thesis, I propose statistically optimal estimators to cancel microphone bleed offline in the mixing and production stage. First, a calibration stage is proposed, where one instrument is played at a time and recorded by all the microphones. This single-input, multiple-output (SIMO) system is used to estimate an approximate relative transfer function matrix, which represents the acoustic path from each source to each microphone and encodes the room response, as well as the mic directivity and source radiation patterns. A convex cost function is derived in the time-frequency domain that simultaneously optimizes the sources and the relative transfer function matrix, which is assumed to be time-invariant. It is shown that minimizing this cost function gives the Maximum Likelihood (ML) estimate when the microphone signals are assumed to be normally distributed. The ML estimator is extended to include a priori statistics of the sources, and the Maximum Aposteriori
Probability (MAP) estimator is derived. The proposed methods are tested on a dataset of anechoic string quartet recordings in a shoebox room, and a drum kit recorded in the studio. Evaluation is done using the PEASS toolbox and a listening test.
</p>


<h3>Publications</h3>
<ul>
	<li><i><a href = "https://www.aes.org/e-lib/browse.cfm?elib=21064">Microphone cross-talk cancellation in ensemble recordings with Maximum Likelihood Estimation</a></i> - O. Das et al. in Proc. of 150th Audio Engineering Society Convention, AES 2021.</li>
	<li><i><a href = "https://ccrma.stanford.edu/~orchi/Documents/odas_thesis_final.pdf">Close-microphone cross-talk cancellation in ensemble recordings with statistical estimation</a></i> - O. Das, PhD thesis, Stanford University, 2021.</li>
</ul>