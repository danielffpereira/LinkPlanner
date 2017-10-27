#include <algorithm>
#include <complex>

#include "netxpto.h"
#include "discarder.h"


void Discarder::initialize(void){

	outputSignals[0]->symbolPeriod = inputSignals[0]->symbolPeriod;;
	outputSignals[0]->samplingPeriod = inputSignals[0]->samplingPeriod;
}


bool Discarder::runBlock(void){
	int ready = inputSignals[0]->ready();
	int space = outputSignals[0]->space();
	int process = min(ready, space);

	if (process == 0) return false;

	for (int i = 0; i < ready; i++) {
		t_binary in1;
		inputSignals[1]->bufferGet(&in1);
		if ( in1 == 0)
		{
			outputSignals[0]->bufferPut(0);
		}
		else if ( in1 == 1)
		{
			outputSignals[0]->bufferPut(1);
		}
		else if ( in1 == 0)
		{
			outputSignals[0]->bufferPut(0);
		}
		else if ( in1 == 1)
		{
			outputSignals[0]->bufferPut(1);
		}
	}

	return true;
}
