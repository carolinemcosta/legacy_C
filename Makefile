all: scar patients vols tetvols classify add-fec map-scar-uvc
	
scar:
	g++ -o scarStudyPlsGrad scarStudy-pls-grad.c
patients:
	g++ -g -o patients-final post-process-patients-final.c -lm
vols:
	g++ -g -o vols computeVolumeRatios.c -lm
tetvols:
	g++ -g -o tetvols computeTetVolumes.c -lm
classify:
	g++ -g -o classify classify-stim.c -lm
add-fec:
	g++ -g -o add-fec add-fec.c -lm	
map-scar-uvc:
	g++ -g -o map-scar-uvc map-scar-to-epi-uvc.c -lm
