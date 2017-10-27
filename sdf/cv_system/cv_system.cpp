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
# include "optical_hybrid.h"
# include "ideal_amplifier.h"
# include "add.h"
# include "white_noise.h" 

int main() {

	// #####################################################################################################
	// #################################### System Input Parameters ########################################
	// #####################################################################################################

	int numberOfBitsReceived(-1);
	int numberOfBitsGenerated(30000);
	int samplesPerSymbol(16);
	int pLength = 0;
	double bitPeriod = 20e-12;
	double rollOffFactor = 0.3;
	double d = 15;
	double T = pow(10,-0.02*d);

	//double t = (T + 2.289) / 3.291; // Xi=0.01

	double t = (T + .2444) / 1.2458; // Xi=0.10
	//double t = (T + .2227) / 1.2240; // Xi=0.11
	//double t = (T + .2046) / 1.2059; // Xi=0.12


	double signalOutputPower = T * 1.0252e-7 / 2; // Half-a-photon energy
	double localOscillatorPower_dBm = 0;
	double localOscillatorPhase = 0;
	vector<t_iqValues> iqAmplitudeValues = { { 1.0, -1.0 },{ -1.0, 1.0 },{ -1.0, -1.0 },{ 1.0, 1.0 } };
	double responsivity = 1;
	double amplification = 1e6;
	double SNU = 102.8412;
	double xi = .10;
	double x = ( xi - ( -0.02443 ) ) / 7.95;
	double electricalNoiseAmplitude = SNU * x * t;
	int samplesToSkip = 8 * samplesPerSymbol;
	int bufferLength = 512;
	bool shotNoise(true);

	// #####################################################################################################
	// ########################### Signals Declaration and Inicialization ##################################
	// #####################################################################################################

	OpticalSignal S1("S1.sgn");
	S1.setBufferLength(bufferLength);

	OpticalSignal S2("S2.sgn");
	S2.setBufferLength(bufferLength);

	OpticalSignal S3("S3.sgn");
	S3.setBufferLength(bufferLength);

	OpticalSignal S4("S4.sgn");
	S4.setBufferLength(bufferLength);

	OpticalSignal S5("S5.sgn");
	S5.setBufferLength(bufferLength);

	OpticalSignal S6("S6.sgn");
	S6.setBufferLength(bufferLength);

	TimeContinuousAmplitudeContinuousReal S7("S7.sgn");
	S7.setBufferLength(bufferLength);

	TimeContinuousAmplitudeContinuousReal S71("S71.sgn");
	S71.setBufferLength(bufferLength);

	TimeContinuousAmplitudeContinuousReal S72("S72.sgn");
	S72.setBufferLength(bufferLength);

	TimeContinuousAmplitudeContinuousReal S8("S8.sgn");
	S8.setBufferLength(bufferLength);

	TimeDiscreteAmplitudeContinuousReal S9("S9.sgn");
	S9.setBufferLength(bufferLength);

	TimeContinuousAmplitudeContinuousReal S10("S10.sgn");
	S10.setBufferLength(bufferLength);

	TimeContinuousAmplitudeContinuousReal S101("S101.sgn");
	S101.setBufferLength(bufferLength);

	TimeContinuousAmplitudeContinuousReal S102("S102.sgn");
	S102.setBufferLength(bufferLength);

	TimeContinuousAmplitudeContinuousReal S11("S11.sgn");
	S11.setBufferLength(bufferLength);

	TimeDiscreteAmplitudeContinuousReal S12("S12.sgn");
	S12.setBufferLength(bufferLength);

	// #####################################################################################################
	// ########################### Blocks Declaration and Inicialization ###################################
	// #####################################################################################################

	MQamTransmitter B1{ vector<Signal*> {}, vector<Signal*> {&S1} };
	B1.setNumberOfBits(numberOfBitsGenerated);
	B1.setOutputOpticalPower(signalOutputPower);
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

	OpticalHybrid B3{ vector<Signal*> {&S1, &S2}, vector<Signal*> {&S3, &S4, &S5, &S6} };
	
	Photodiode B4{ vector<Signal*> {&S3, &S4}, vector<Signal*> {&S7} };
	B4.useNoise(shotNoise);
	B4.setResponsivity(responsivity);
	B4.generatorAmp1.seed(3585632034);
	B4.generatorAmp2.seed(9873456547);

	TI_Amplifier B5{ vector<Signal*> {&S7}, vector<Signal*> {&S8} };
	B5.setGain(amplification);
	B5.setElectricalNoiseSpectralDensity(electricalNoiseAmplitude);
	B5.setNoiseSeed(2228957057);

	Sampler B6{ vector<Signal*> {&S8}, vector<Signal*> {&S9} };
	B6.setSamplesToSkip(samplesToSkip);
	
	Photodiode B7{ vector<Signal*> {&S5, &S6}, vector<Signal*> {&S10} };
	B7.useNoise(shotNoise);
	B7.setResponsivity(responsivity);
	B7.generatorAmp1.seed(3337066948);
	B7.generatorAmp2.seed(6843974536);

	TI_Amplifier B8{ vector<Signal*> {&S10}, vector<Signal*> {&S11} };
	B8.setGain(amplification);
	B8.setElectricalNoiseSpectralDensity(electricalNoiseAmplitude);
	B8.setNoiseSeed(3435973836);

	Sampler B9{ vector<Signal*> {&S11}, vector<Signal*> {&S12} };
	B9.setSamplesToSkip(samplesToSkip);

	Sink B10{ vector<Signal*> { &S9 }, vector<Signal*> {} };
	B10.setNumberOfSamples(samplesPerSymbol*numberOfBitsGenerated);
	B10.setDisplayNumberOfSamples(true);

	Sink B11{ vector<Signal*> { &S12 }, vector<Signal*> {} };
	B11.setNumberOfSamples(samplesPerSymbol*numberOfBitsGenerated);
	B11.setDisplayNumberOfSamples(true);


	// #####################################################################################################
	// ########################### System Declaration and Inicialization ###################################
	// #####################################################################################################

	System MainSystem{ vector<Block*> { &B1, &B2, &B3, &B4, &B5, &B6, &B7, &B8, &B9, &B10, &B11} };

	// #####################################################################################################
	// #################################### System Run #####################################################
	// #####################################################################################################

	MainSystem.run();

	return 0;

}
