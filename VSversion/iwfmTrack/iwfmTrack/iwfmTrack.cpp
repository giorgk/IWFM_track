// iwfmTrack.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <functional>
#include <iostream>
#include <chrono>
#include <thread>

#include "iwfmOptions.h"
#include "iwfmReader.h"
#include "iwfmTree.h"
#include "wells.h"
#include "ParticleTack.h"

int main(int argc, char *argv[])
{
	iwfm_track_options iwfm_opt;
	bool tf = readInputParameters(argc, argv, iwfm_opt);
	if (!tf)
		return 0;

	VelocityField HF(iwfm_opt.NRMLfile, iwfm_opt.FACEZfile, iwfm_opt.HFLOWfile, iwfm_opt.radius, iwfm_opt.Nsearch, iwfm_opt.Nminp, iwfm_opt.threshold);
	VelocityField VF(iwfm_opt.BCfile, iwfm_opt.BCZfile, iwfm_opt.VFLOWfile, iwfm_opt.radius, iwfm_opt.Nsearch, iwfm_opt.Nminp, iwfm_opt.threshold);
	iwfmMesh MSH(iwfm_opt.XYfile, iwfm_opt.MSHfile, iwfm_opt.STRATfile, iwfm_opt.threshold);

	inputWells wells(iwfm_opt);
	wells.distributeParticles(MSH);

	pTrack pt(HF, VF, MSH, iwfm_opt);

	auto start = std::chrono::high_resolution_clock::now();

	bool bThreaded = iwfm_opt.nThreads > 0;
	if (bThreaded) {
		std::vector<std::thread> T;
		for (int i = 0; i < iwfm_opt.nThreads; ++i) {
			T.push_back(std::thread(&pTrack::trace_with_threads, std::ref(pt), i, std::ref(wells.initStreamlines)));
		}
		for (int i = 0; i < iwfm_opt.nThreads; ++i) {
			T[i].join();
		}
	}
	else {
		pt.trace(wells.initStreamlines);
	}
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Particle tracking finished in " << elapsed.count() << std::endl;

	wells.WriteStreamlines();

	return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
