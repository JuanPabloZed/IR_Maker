-- IR Maker --
Developped by Nathan Zwahlen, Benjamin Quiédeville & Hugo Perrier.

To use the software, extract the "IR Maker.zip" file (takes roughly 250MB when extracted) and
run the "IR Maker.exe" file.

This application is made to create impulse responses (IRs) using the Exponential Sine Sweep (ESS)
method (for more info, check Angelo Farina's ESS method). This method consists of the following : 

	1- Generate an ESS. If you don't have one, you can create one with IR Maker : 
	just press the Sweep Generator button and follow the instructions.
	
	2- Play the ESS through your signal chain (speaker if you want a room/reverb IR, amp
	and/or cab if you want a cab IR, or any effects chain, hardware or in a DAW).
	
	3- Record the full response of the signal chain to the ESS. Make sure not to cut any
	data by ending too soon the recording.
	
	4- Export the response in WAV format (be sure to export in MONO if it's a cab IR).
	
	5- Run IR Maker, follow the instructions in the software.
	
	6- Drop the IR in your DAW if you want to edit it for an even cleaner IR.

And voilà ! Here is your custom-made IR, usable in any DAW, convolution reverb/ampsim/IR loader
plugin.

/!\ WARNING - Sample rate matching /!\ 
Be sure to use the same sample rate for the ESS, the response and the IR to get the cleanest IR.
Otherwise, the IR might sound unexpectedly weird. But if you want to experiment, that's a nice
way to start !

Thank you for using our software !