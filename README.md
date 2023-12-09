-- IR Maker --
Developped by Nathan Zwahlen, Benjamin Quiédeville & Hugo Perrier.

This application is made to create impulse responses (IRs) using the Exponential Sine Sweep (ESS)
method (for more info, check Angelo Farina's ESS method). This method consists of the following : 

	1- Generate an ESS. If you don't have one, you can create one with IR Maker.
	
	2- Play the ESS through your signal chain (speaker if you want a room/reverb IR, amp
	and/or cab if you want a cab IR).
	
	3- Record the full response of your room or amp/cab to the ESS. Make sure not to cut any
	data by ending too soon the recording.
	
	4- Export the response in WAV format (be sure to export in mono if it's a cab IR).
	
	5- Run IR Maker, follow the instructions in the UI.
	
	6- Drop the IR in your DAW if you want to edit it for an even cleaner IR.

And voilà ! Here is your custom-made IR, usable in any DAW, convolution reverb/ampsim/IR loader
plugin.

/!\ WARNING /!\ 
Be sure to use the same sample rate for the ESS, the response and the IR to get the cleanest IR.

/!\ COMPILING INTO .exe /!\

If you want to make a .exe out of this code, you can by installing the different necessary modules contained in requirements.txt, and then typing the following command in cmd: 
pyinstaller --noconfirm --onefile --windowed "path\of\the\output.exe\file"


Thank you for using our application !
