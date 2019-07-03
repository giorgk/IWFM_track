#pragma once
#include <vector>
#include "iwfmTree.h"
#include "iwfmOptions.h"

struct vec3D {
	double x;
	double y;
	double z;

	vec3D(double x = 0, double y = 0, double z = 0)
		:
		x(x), y(y), z(z)
	{}
	
	vec3D operator*(const vec3D& a) const {
		return vec3D(a.x*x, a.y*y, a.z*z);
	}

	vec3D operator*(const double& a) const {
		return vec3D(a*x, a*y, a*z);
	}

	vec3D operator+(const vec3D& a) const {
		return vec3D(a.x+x, a.y+y, a.z+z);
	}

	vec3D operator-(const vec3D& a) const {
		return vec3D(a.x - x, a.y - y, a.z - z);
	}

	double len() {
		return std::sqrt(x*x + y * y + z * z);
	}
};

struct particle {
	vec3D p;
	vec3D v;
	double age = 0;
	bool hasVelocity = false;
};

struct streamline {
	int Eid;
	int Sid;
	std::vector<particle> SL;
	int exitflag = 0;
	void calcAge() {
		if (SL.size() > 1) {
			double v = SL[SL.size() - 2].v.len();
			double a = SL[SL.size() - 2].age;
			double d = (SL[SL.size() - 1].p - SL[SL.size() - 2].p).len();
			SL[SL.size() - 1].age = a + d / v;
		}
	}
};


class pTrack {
public:
	pTrack(VelocityField&	Hin,
		VelocityField&		Vin,
		iwfmMesh&			MSHin,
		iwfm_track_options&	optin);
	void trace(std::vector<streamline>& Streamlines);

private:
	VelocityField		H;
	VelocityField		V;
	iwfmMesh			MSH;
	iwfm_track_options	opt;
	double					initstep;

	void traceStreamline(streamline& S);

	void calcVel(particle &p);
	int isinDomain(particle p);
	particle findnextpoint(particle p);
};

pTrack::pTrack(VelocityField&	Hin,
	VelocityField&		Vin,
	iwfmMesh&			MSHin,
	iwfm_track_options&	optin)
	:
	H(Hin),
	V(Vin),
	MSH(MSHin),
	opt(optin) {
	initstep = opt.stepSize;
}

void pTrack::trace(std::vector<streamline>& Streamlines) {
	std::vector<streamline>::iterator its;
	for (its = Streamlines.begin(); its != Streamlines.end(); ++its) {
		traceStreamline(*its);
	}
}

void pTrack::traceStreamline(streamline& S) {
	S.SL.at(0);
	opt.stepSize = initstep;
	std::cout << "Tracing Eid: " << S.Eid << ", Sid: " << S.Sid << std::endl;
	int exitflag = isinDomain(S.SL.back());
	int count_iter = 0;
	while (exitflag == 0) {
		if (!S.SL.back().hasVelocity) {
			calcVel(S.SL.back());
		}
		S.SL.push_back(findnextpoint(S.SL.back()));
		S.calcAge();
		exitflag = isinDomain(S.SL.back());

		count_iter++;
		if (count_iter > opt.maxSteps) {
			exitflag = -8;
			break;
		}
	}
	S.exitflag = exitflag;
}

int pTrack::isinDomain(particle p) {
	std::vector<double> elev;
	if (!MSH.getPntStrat(p.p.x, p.p.y, elev)) {
		return 2;
	}
	else {
		if (p.p.z > elev.front()) {
			return 1;
		}
		
		else if (p.p.z < elev.back()) {
			return -9;
		}
		else {
			return 0;
		}
	}
}

void pTrack::calcVel(particle &p) {
	std::vector<double> vxy;
	H.XYcalc(p.p.x, p.p.y, p.p.z, vxy);
	double vz = V.Zcalc(p.p.x, p.p.y, p.p.z);
	p.v.x = opt.dir * vxy[0];
	p.v.y = opt.dir * vxy[1];
	p.v.z = opt.dir * vz;
	p.hasVelocity = true;
}

particle pTrack::findnextpoint(particle p) {
	vec3D k1, k2, k3, k4, k5, k6;
	k1 = p.v;
	double vmag = k1.len();
	double m = opt.stepSize / vmag;

	particle p2, p3, p4, p5, p6, yn, zn;

	p2.p= p.p + k1 * m*(1.0 / 4.0);
	calcVel(p2);
	k2 = p2.v;

	p3.p = p.p + ( k1*(3.0 / 32.0) + k2 * (9.0 / 32.0) ) * (3.0 / 8.0)*m;
	calcVel(p3);
	k3 = p3.v;

	p4.p = p.p + ( k1*(1932.0 / 2197.0) + k2*(-(7200.0 / 2197.0)) + k3* (7296.0 / 2197.0) ) * (12.0 / 13.0)*m;
	calcVel(p4);
	k4 = p4.v;

	p5.p = p.p + ( k1*(439.0 / 216.0) + k2*(-8.0) + k3*(3680.0 / 513.0) + k4*(-(845.0 / 4104.0)) ) * (1.0)*m;
	calcVel(p5);
	k5 = p5.v;

	p6.p = p.p + ( k1*(-(8.0 / 27.0)) + k2*2.0 + k3*(-(3544.0 / 2565.0)) + k4* (1859.0 / 4104.0) + k5*(-(11.0 / 40.0)) ) * (1.0 / 2.0)*m;
	calcVel(p6);
	k6 = p6.v;

	yn.p = p.p + ( k1*(25.0 / 216.0) + k3*(1408.0 / 2565.0) + k4* (2197.0 / 4101.0) + k5*(-(1.0 / 5.0)) ) *m;
	zn.p = p.p + ( k1*(16.0 / 135.0) + k3*(6656.0 / 12825.0) + k4*(28561.0 / 56430.0) + k5*(-(9.0 / 50.0)) + k6* (2.0 / 55.0) ) *m;

	double R = (zn.p - yn.p).len();
	double delta = 0.84*pow(opt.tolStepSize / R, 0.25);
	opt.stepSize = delta * opt.tolStepSize*opt.stepSize;
	if (R < opt.tolStepSize) {
		
		if (opt.stepSize < opt.minStepSize)
			opt.stepSize = opt.minStepSize;
		if (opt.stepSize > opt.maxStepSize)
			opt.stepSize = opt.maxStepSize;
		return zn;
	}
	else {
		return findnextpoint(p);
	}
	particle pnext;
	return pnext;
}
