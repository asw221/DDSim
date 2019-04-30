
DDSim.exe: DDSim.cpp DD-Arrangement.h DD-Candidate.h DD-Data.h DD-Donor.h DD-Enums-Functions.h DD-Match.h DD-MatchRun.h DD-Node.h DD-Parameters.h DD-Record.h DD-RNG.h DD-Simulation.h
	clang++ -Wall -g -O0 -std=c++17 -pthread -L$(GUROBI_HOME)/lib -I$(GUROBI_HOME)/include -lgurobi_c++ -lgurobi80 -o DDSim.exe DDSim.cpp


clean:
	rm -rf *.exe




