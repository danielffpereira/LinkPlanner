# include "netxpto.h"

# include "m_qam_transmitter.h"
# include "i_homodyne_receiver.h"
# include "bit_error_rate.h"
# include "local_oscillator.h"
# include "balanced_beam_splitter.h"
# include "photodiode_old.h"
# include "ti_amplifier.h"
# include "sampler.h"
# include "sink.h"

int main(){

	// #####################################################################################################
	// #################################### System Input Parameters ########################################
	// #####################################################################################################

	int numberOfBitsReceived(-1);
	int numberOfBitsGenerated(300);
	int samplesPerSymbol(16);
	int pLength = 5;
	double bitPeriod = 20e-12;
	double rollOffFactor = 0.3;
	double signalOutputPower_dBm = -20;
	double localOscillatorPower_dBm = -20;
	double localOscillatorPhase = 0;
	vector<t_iqValues> iqAmplitudeValues = { {-1, 0} , {1, 0} };
	array<t_complex, 4> transferMatrix = { { 1 / sqrt(2), 1 / sqrt(2), 1 / sqrt(2), -1 / sqrt(2)} };
	double responsivity = 1;
	double amplification = 1;
	double electricalNoiseAmplitude = 0;// 15.397586549153788;
	int samplesToSkip = 8 * samplesPerSymbol;
	int bufferLength = 512;
	bool shotNoise(false);
		
	// #####################################################################################################
	// ########################### Signals Declaration and Inicialization ##################################
	// #####################################################################################################

	Binary S0("S0.sgn");
	S0.setBufferLength(bufferLength);

	OpticalSignal S1("S1.sgn");
	S1.setBufferLength(bufferLength);

	OpticalSignal S2("S2.sgn");
	S2.setBufferLength(bufferLength);

	OpticalSignal S3("S3.sgn");
	S3.setBufferLength(bufferLength);

	OpticalSignal S4("S4.sgn");
	S4.setBufferLength(bufferLength);

	TimeContinuousAmplitudeContinuousReal S5("S5.sgn");
	S5.setBufferLength(bufferLength);

	TimeContinuousAmplitudeContinuousReal S6("S6.sgn");
	S6.setBufferLength(bufferLength);

	TimeDiscreteAmplitudeContinuousReal S7("S7.sgn");
	S7.setBufferLength(bufferLength);

	Binary S8("S8.sgn");
	S8.setBufferLength(bufferLength);

	Binary S9("S9.sgn");
	S9.setBufferLength(bufferLength);

	// #####################################################################################################
	// ########################### Blocks Declaration and Inicialization ###################################
	// #####################################################################################################

	MQamTransmitter B1{ vector<Signal*> {}, vector<Signal*> {&S1, &S0} };
	B1.setNumberOfBits(numberOfBitsGenerated);
	B1.setOutputOpticalPower_dBm(signalOutputPower_dBm);
	B1.setMode(PseudoRandom);
	B1.setBitPeriod(bitPeriod);
	B1.setPatternLength(pLength);
	B1.setIqAmplitudes(iqAmplitudeValues);
	B1.setNumberOfSamplesPerSymbol(samplesPerSymbol);
	B1.setRollOffFactor(rollOffFactor);
	B1.setSaveInternalSignals(true);
	B1.setPulseShaperFilter(RaisedCosine);
	B1.setSeeBeginningOfImpulseResponse(false);



	LocalOscillator B2{ vector<Signal*> { }, vector<Signal*> { &S2 } };
	B2.setOpticalPower_dBm(localOscillatorPower_dBm);
	B2.setPhase(localOscillatorPhase);
	B2.setSamplingPeriod(bitPeriod / samplesPerSymbol);
	B2.setSymbolPeriod(bitPeriod);

	BalancedBeamSplitter B3{ vector<Signal*> {&S1, &S2}, vector<Signal*> {&S3, &S4} };
	B3.setTransferMatrix(transferMatrix);

	Photodiode B4{ vector<Signal*> {&S3, &S4}, vector<Signal*> {&S5} };
	B4.useNoise(shotNoise);
	B4.setResponsivity(responsivity);

	TI_Amplifier B5{ vector<Signal*> {&S5}, vector<Signal*> {&S6} };
	B5.setGain(amplification);
	B5.setElectricalNoiseSpectralDensity(electricalNoiseAmplitude);
	B5.setSaveInternalSignals(true);
	B5.setSeeBeginningOfImpulseResponse(false);
	B5.setImpulseResponseLength(16);
	B5.setRollOffFactor(0);
	B5.usePassiveFilterMode(true);

	Sampler B6{ vector<Signal*> {&S6}, vector<Signal*> {&S7} };
	B6.setSamplesToSkip(samplesToSkip);
	
	BitDecider B7{ vector<Signal*> {&S7}, vector<Signal*> {&S8} };

	BitErrorRate B8{ vector<Signal*> {&S0, &S8}, vector<Signal*> {&S9} };

	Sink B9{ vector<Signal*> { &S9 }, vector<Signal*> {} };
	B9.setNumberOfSamples(samplesPerSymbol*numberOfBitsGenerated);
	B9.setDisplayNumberOfSamples(true);


	// #####################################################################################################
	// ########################### System Declaration and Inicialization ###################################
	// #####################################################################################################

	System MainSystem{ vector<Block*> { &B1, &B2, &B3, &B4, &B5, &B6, &B7, &B8, &B9} };

	// #####################################################################################################
	// #################################### System Run #####################################################
	// #####################################################################################################

	MainSystem.run();

	return 0;

}
