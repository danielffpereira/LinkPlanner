#include <algorithm>
#include <complex>
#include <iostream>
#include <fstream>

#include "netxpto.h"
#include "sampler.h"


void Sampler::initialize(void){

	firstTime = false;

	outputSignals[0]->setSymbolPeriod(inputSignals[0]->getSymbolPeriod());
	outputSignals[0]->setSamplingPeriod(inputSignals[0]->getSymbolPeriod());
}


bool Sampler::runBlock(void){

	ofstream myfile2;
	myfile2.open("translate.txt", fstream::app);

	int ready = inputSignals[0]->ready();

	int samplesPerSymbol = inputSignals[0]->getSamplesPerSymbol();
	
	if (firstPass || repeatedPass) 
	{

		if (firstPass)
		{
			samplesToSkip = samplesPerSymbol * 8 + 336;
			aux1 = true;
		}

		firstPass = false;

		int process = min(ready, samplesToSkip);


		for (int k = 0; k < process; k++) {
			t_real in;
			inputSignals[0]->bufferGet(&in);
		}

		samplesToSkip = samplesToSkip - process;
		repeatedPass = false;
		if (samplesToSkip != 0)
		{
			repeatedPass = true;
		}
		ready = inputSignals[0]->ready();
	}

	int space = outputSignals[0]->space();
	int process = min(ready, space);
	
	
	if (process == 0){
		myfile2.close();
		return false;
	}
	if (samplesToSkip == 0)
	{
		for (int k = 0; k < process; k++) {
			t_real in;
			inputSignals[0]->bufferGet(&in);
			if (count % samplesPerSymbol == 0) {

<<<<<<< HEAD
			inputSignals[1]->bufferGet(&inClock);
			inputSignals[0]->bufferGet(&inSignal);

			if (inClock == 1.0) {

				inSignal = inSignal; // (.5*sqrt(outputOpticalPower)); // to normalize the signal to 1
				outputSignals[0]->bufferPut(inSignal);

=======
				outputSignals[0]->bufferPut((t_real)in);
>>>>>>> develop
			}
			myfile2 << in << "\n";
		}

	}
<<<<<<< HEAD

	return 0;

};
=======
	return true;
}
>>>>>>> develop
