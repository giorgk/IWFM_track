#pragma once

#include <iostream>
#include <string>
#include <fstream>

#include "iwfmOptions.h"
#include "ParticleTack.h"
#include "iwfmTree.h"

struct well {
	double x;
	double y;
	double d; // depth from ground surface to water table
	double s; // screen length
};

class inputWells {
public:
	inputWells(iwfm_track_options & opt);
	void distributeParticles(iwfmMesh & MSH);

	std::vector<streamline> initStreamlines;
	void WriteStreamlines();
private:
	void readwellFile(std::string filename);
	void readParticleFile(std::string filename);
	unsigned int Nwells;
	std::vector<well> Wells;
	int Npart;
	double radius;
	std::string outfile;
	int pntPrec;
	int velPrec;
	void writeEndpnts();
	bool printendpnts;
};


inputWells::inputWells(iwfm_track_options & opt) {
	readwellFile(opt.WELLSfile);
	readParticleFile(opt.PARTICLEfile);
	Npart = opt.Npart;
	radius = opt.wellRadius;
	outfile = opt.OUTfile;
	pntPrec = opt.pntPrecision;
	velPrec = opt.velPrecision;
	printendpnts = opt.bEndpnts;
}

void inputWells::readwellFile(std::string filename) {
	if (filename.empty())
		return;

	std::ifstream datafile(filename.c_str());
	if (!datafile.good()) {
		std::cout << "Can't open the file" << filename << std::endl;
	}
	else {
		char buffer[1024];
		datafile.getline(buffer, 1024);
		std::istringstream inp(buffer);
		inp >> Nwells;
		for (unsigned int i = 0; i < Nwells; ++i) {
			datafile.getline(buffer, 1024);
			std::istringstream inp(buffer);
			well w;
			inp >> w.x;
			inp >> w.y;
			inp >> w.d;
			inp >> w.s;

			Wells.push_back(w);
		}
	}
	datafile.close();
}

void inputWells::readParticleFile(std::string filename) {
	if (filename.empty())
		return;

	std::ifstream datafile(filename.c_str());
	if (!datafile.good()) {
		std::cout << "Can't open the file" << filename << std::endl;
	}
	else {
		int Nparticles;
		double x, y, z;
		int Eid, Sid;
		char buffer[1024];
		datafile.getline(buffer, 1024);
		std::istringstream inp(buffer);
		inp >> Nparticles;
		for (unsigned int i = 0; i < Nparticles; ++i) {
			datafile.getline(buffer, 1024);
			std::istringstream inp(buffer);
			inp >> Eid;
			inp >> Sid;
			inp >> x;
			inp >> y;
			inp >> z;
			streamline s;
			s.Eid = Eid;
			s.Sid = Sid;
			vec3D pos;
			pos.x = x;
			pos.y = y;
			pos.z = z;
			particle p;
			p.p = pos;
			s.SL.push_back(p);
			initStreamlines.push_back(s);
		}
	}
}

void inputWells::distributeParticles(iwfmMesh & MSH) {

	std::vector<well>::iterator itw;

	int cnt_wells = 1;
	for (itw = Wells.begin(); itw != Wells.end(); ++itw ) {
		std::vector<double> stratpnt;
		if (MSH.getPntStrat(itw->x, itw->y, stratpnt)) {
			double gse = stratpnt.at(0);
			double base = stratpnt.at(stratpnt.size() - 1);
			double wellTop = gse - itw->d;
			double wellBot = wellTop - itw->s;
			for (int i = 0; i < Npart; ++i) {
				double t = static_cast<double>(i + 1);
				particle p;
				vec3D pos;
				pos.x = radius * std::cos(t) + itw->x;
				pos.y = radius * std::sin(t) + itw->y;
				pos.z = (t * (1 / static_cast<double>(Npart)))*(wellTop - wellBot) + wellBot;
				p.p = pos;
				//std::cout << std::setprecision(5) << std::fixed << p.x << " " << p.y << " " << p.z << std::endl;

				streamline s;
				s.Eid = cnt_wells;
				s.Sid = i;
				s.SL.push_back(p);
				initStreamlines.push_back(s);
			}
		}
		++cnt_wells;
	}
}

void inputWells::WriteStreamlines() {
	if (printendpnts) {
		writeEndpnts();
	}
	else {
		std::cout << "Writing streamlines" << std::endl;
		std::ofstream outstream;
		outstream.open(outfile.c_str());

		outstream << initStreamlines.size() << std::endl;
		std::vector<streamline>::iterator its;
		for (its = initStreamlines.begin(); its != initStreamlines.end(); ++its) {
			outstream << its->Eid << " " << its->Sid << " " << its->exitflag << " " << its->SL.size() << " " << std::endl;
			for (unsigned int i = 0; i < its->SL.size(); ++i) {
				outstream << std::setprecision(pntPrec) << std::fixed
					<< its->SL[i].p.x << " " 
					<< its->SL[i].p.y << " " 
					<< its->SL[i].p.z << " " 
					<< std::setprecision(velPrec) << std::fixed
					<< its->SL[i].v.x << " " 
					<< its->SL[i].v.y << " " 
					<< its->SL[i].v.z << " " 
					<< its->SL[i].age << std::endl;
			}
		}
		outstream.close();
	}
}

void inputWells::writeEndpnts() {
	std::cout << "Writing Endpoints of streamlines" << std::endl;
	std::ofstream outstream;
	outstream.open(outfile.c_str());

	outstream << initStreamlines.size() << std::endl;
	std::vector<streamline>::iterator its;
	for (its = initStreamlines.begin(); its != initStreamlines.end(); ++its) {
		outstream << its->Eid << " " << its->Sid << " " << its->exitflag << " " 
			<< std::setprecision(pntPrec) << std::fixed
			<< its->SL.back().p.x << " " 
			<< its->SL.back().p.y << " " 
			<< its->SL.back().p.z << " " 
			<< its->SL.back().age << " " << std::endl;
	}
	outstream.close();
}
