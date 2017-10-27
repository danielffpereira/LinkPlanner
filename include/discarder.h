# ifndef PROGRAM_INCLUDE_DISCARDER_H_
# define PROGRAM_INCLUDE_DISCARDER_H_

# include "netxpto.h"
# include <vector>

// Adjusts the starting point of the discretized signal
class Discarder : public Block {
public:
	Discarder(vector<Signal *> &InputSig, vector<Signal *> &OutputSig) :Block(InputSig, OutputSig) {};
	
	void initialize(void);
	bool runBlock(void);

private:

	int aux = 0;

};


#endif // !PROGRAM_INCLUDE_DISCARDER_H_
