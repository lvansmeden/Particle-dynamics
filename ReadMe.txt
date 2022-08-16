This folder contains three executable MATLAB-codes which can be used to generate motion patterns of the f-BPM System:

- f_BPM_OnlyFreeDiffusion
	This code only simulates free diffusive motion of a particle in the f-BPM system.
	No protrusions are put on the sphere's surface (= smooth).

- f_BPM_MotionSimulator
	This code simulates motion patterns that include switching between unbound, single-bound and double-bound motion
	First section of the code serves as input and can be altered
	
	When code asks to reducte the integration time to resolve stability issues: use lines 143-145 of the script
	Output is saved in an Excel file
	In Excel file, the x-position (m), y-position (m) and number of bonds are saved per frame.

- f_BPM_MotionSimulator_Bunch
	This script basically runs the same script as the script above, but enables the simulation of multiple particles after another
	Allows for including variability in protrusions on particle surface or distance between binders on the substrate needs.
	Output is again saved in an Exel file.
	Output looks like
		[1	2	3	4	5	6 	..... 		]
		[x-position particle 1	y-position particle 1	state particle  1	x-position particle 2	y-position particle 2	state particle  2 ..... ]
		etc.