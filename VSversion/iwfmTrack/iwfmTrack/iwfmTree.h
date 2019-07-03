#pragma once

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/foreach.hpp>

#include "iwfmReader.h"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<float, 2, bg::cs::cartesian> point_b;
typedef bg::model::polygon<point_b> poly_b;
typedef std::pair<point_b, unsigned int> value;

class boostTree {
public:
	boostTree();
	void buildTree(iwfmMatrix2D<double>& XY);
	//bool searchTree(float x, float y, std::vector<std::pair<double, int>> &id);
	bool searchTree(float x, float y, std::vector<int> &id);
	void setNminPnts(int depth);
	void setDistance(double d);
	void setNsearch(int N);

private:
	bgi::rtree <value, bgi::quadratic<16> > rTree;
	int NminPnts;
	int Nsearch;
	double dist;
};

boostTree::boostTree(){}

void boostTree::buildTree(iwfmMatrix2D<double>& XY) {
	std::cout << "Building tree..." << std::endl;
	for (int i = 0; i < XY.nrows(); ++i) {
		rTree.insert(std::make_pair(point_b(
			static_cast<float>(XY.get(i, 0)),
			static_cast<float>(XY.get(i, 1))
		), i));
	}
	std::cout << "Building tree finish" << std::endl;
}

void boostTree::setNsearch(int N) {
	Nsearch = N;
}

void boostTree::setNminPnts(int N) {
	NminPnts = N;
}

void boostTree::setDistance(double d) {
	dist = d;
}

bool boostTree::searchTree(float x, float y, std::vector<int> &id) {
	id.clear();
	std::vector<value> result_n;
	rTree.query(bgi::nearest(point_b(x, y), Nsearch), std::back_inserter(result_n));

	std::vector<std::pair<double, int>> sortedRes;
	BOOST_FOREACH(value const& v, result_n) {
		float xn = v.first.get<0>();
		float yn = v.first.get<1>();
		float dst = std::sqrt((x - xn)*(x - xn) + (y - yn) * (y - yn));
		sortedRes.push_back(std::make_pair(dst, v.second));
		//std::cout << bg::wkt<point_b>(v.first) << " - " << v.second << " - " << dst << std::endl;
	}

	std::sort(sortedRes.begin(), sortedRes.end());
	for (unsigned int i = 0; i < sortedRes.size(); ++i) {
		if (sortedRes[i].first < dist) {
			id.push_back(sortedRes[i].second);
		}
		else {
			if (id.size() < static_cast<unsigned int>(NminPnts))
				id.push_back(sortedRes[i].second);
			else
				break;
		}
	}

	return false;
}

class VelocityField {
public:
	VelocityField(std::string xynrFile, std::string Zfile, std::string Vfile, double rad, int Nsearch, int NminPnts, double tol);
	void XYcalc(double x, double y, double z, std::vector<double>& vxy);
	double Zcalc(double x, double y, double z);
	

private:
	// matrix that holds xy and normal vectors
	iwfmMatrix2D<double> XYNR;
	// matrix that holds elevation of face centers
	iwfmMatrix2D<double> Z;
	// matrix that holds the velocity magnitude 
	iwfmMatrix2D<double> V;

	double thres;

	boostTree XYtree;

	void findZpoints(std::vector<std::pair<int, int>> &pnts, double z, std::vector<int>& id);

	void IDWinterpolateXY(std::vector<double> &vxy, double x, double y, double z, std::vector<std::pair<int, int>> &pnts);
	double IDWinterpolateZ(double x, double y, double z, std::vector<std::pair<int, int>> &pnts);
};

VelocityField::VelocityField(std::string xynrFile, std::string Zfile, std::string Vfile, double rad, int Nsearch, int NminPnts, double tol) {
	XYNR.readFromFile(xynrFile, false);
	Z.readFromFile(Zfile, false);
	V.readFromFile(Vfile, false);

	XYtree.buildTree(XYNR);
	XYtree.setNminPnts(NminPnts);
	XYtree.setNsearch(Nsearch);
	XYtree.setDistance(rad);
	thres = tol;
}


void VelocityField::XYcalc(double x, double y, double z, std::vector<double>& vxy) {
	if (XYNR.ncols() < 3) {
		std::cerr << "Normals are missing" << std::endl;
		return;
	}
	std::vector<int> ids;
	bool tf = XYtree.searchTree(static_cast<float>(x), static_cast<float>(y), ids);
	std::vector<std::pair<int, int>> pnts;
	findZpoints(pnts, z, ids);
	IDWinterpolateXY(vxy, x, y, z, pnts);
}

double VelocityField::Zcalc(double x, double y, double z) {
	std::vector<int> ids;
	bool tf = XYtree.searchTree(static_cast<float>(x), static_cast<float>(y), ids);
	std::vector<std::pair<int, int>> pnts;
	findZpoints(pnts, z, ids);
	return IDWinterpolateZ(x, y, z, pnts);
}

void VelocityField::findZpoints(std::vector<std::pair<int, int>> &pnts, double z, std::vector<int>& id) {
	pnts.clear();
	for (unsigned int ii = 0; ii < id.size(); ++ii) {
		int i = id[ii];
		for (int j = 0; j < Z.ncols(); ++j) {
			if (j == 0) {
				if (z > Z.get(i,j)) {
					pnts.push_back(std::pair<int, int>(i, j));
					break;
				}
			}
			else {
				if (j == Z.ncols() - 1 && z < Z.get(i, j)) {
					pnts.push_back(std::pair<int, int>(i, j));
					break;
				}
				if (Z.get(i, j-1) > z && Z.get(i, j) < z) {
					pnts.push_back(std::pair<int, int>(i, j - 1));
					pnts.push_back(std::pair<int, int>(i, j));
					break;
				}
			}
		}
	}
}

void VelocityField::IDWinterpolateXY(std::vector<double> &vxy, double x, double y, double z, std::vector<std::pair<int, int>> &pnts) {
	double sumW = 0;
	double sumWV = 0;
	double sumWNX = 0;
	double sumWNY = 0;
	std::vector<std::pair<int, int>>::iterator it;
	for (it = pnts.begin(); it != pnts.end(); ++it) {
		int i = it->first;
		int j = it->second;
		double d = std::sqrt(pow(x - XYNR.get(i,0), 2) + pow(y - XYNR.get(i, 1), 2) + pow(z - Z.get(i,j), 2));
		if (d < thres) {
			vxy.push_back(V.get(i, j) * XYNR.get(i, 2));
			vxy.push_back(V.get(i, j) * XYNR.get(i, 3));
			break;
		}
		else {
			d = 1.0 / d;
			sumW += d;
			sumWV += d * V.get(i, j);
			sumWNX += d * XYNR.get(i, 2);
			sumWNY += d * XYNR.get(i, 3);
		}
	}

	double nx = sumWNX / sumW;
	double ny = sumWNY / sumW;
	double v = sumWV / sumW;
	double L = std::sqrt(nx*nx + ny * ny);
	nx = nx / L;
	ny = ny / L;
	vxy.push_back(v*nx);
	vxy.push_back(v*ny);
}

double VelocityField::IDWinterpolateZ(double x, double y, double z, std::vector<std::pair<int, int>> &pnts) {
	double sumW = 0;
	double sumWV = 0;
	std::vector<std::pair<int, int>>::iterator it;
	for (it = pnts.begin(); it != pnts.end(); ++it) {
		int i = it->first;
		int j = it->second;
		double d = std::sqrt(pow(x - XYNR.get(i, 0), 2) + pow(y - XYNR.get(i, 1), 2) + pow(z - Z.get(i, j), 2));
		double v = 0;
		if (j < V.ncols()) // To account for the fact that the last layer has zero velocity and it is not store in the V matrix.
			v = V.get(i, j);

		if (d < thres) {
			
			return v;
		}
		else {
			d = 1.0 / d;
			sumW += d;
			sumWV += d * v;
		}
	}
	return sumWV / sumW;
}

class iwfmMesh {
public:
	iwfmMesh(std::string XYfile, std::string MSHfile, std::string STRATfile, double tol);
	int findElement(double x, double y);
	bool getPntStrat(double x, double y, std::vector<double> &elev);
	double getXY(int i, int j);
	int getMSH(int i, int j);
	double getZ(int i, int j);
	int Nnodes();
	int Nelem();
	int Nlay();

private:
	iwfmMatrix2D<double> XY;
	iwfmMatrix2D<int> MSH;
	iwfmMatrix2D<double> STRAT;
	void calcBCTree();
	boostTree BCtree;
	double tolerance;
};


iwfmMesh::iwfmMesh(std::string XYfile, std::string MSHfile, std::string STRATfile, double tol) {
	XY.readFromFile(XYfile, false);
	MSH.readFromFile(MSHfile, false);
	STRAT.readFromFile(STRATfile, false);
	tolerance = tol;
	calcBCTree();
}

void iwfmMesh::calcBCTree() {
	std::vector<std::vector<double>> BC;
	double bcx, bcy, n;
	
	for (int i = 0; i < MSH.nrows(); ++i) {
		bcx = 0; bcy = 0;
		n = 4;
		for (int j = 0; j < 4; ++j) {
			if (j == 3 && MSH.get(i, j) == 0) {
				n = 3.0;
				break;
			}
			bcx += XY.get(MSH.get(i, j) - 1, 0);
			bcy += XY.get(MSH.get(i, j) - 1, 1);
			
		}
		bcx = bcx / n;
		bcy = bcy / n;
		BC.push_back(std::vector<double>{bcx, bcy});
	}
	iwfmMatrix2D<double> iwfmBC;
	iwfmBC.fromVector(BC, false);

	BCtree.buildTree(iwfmBC);
	BCtree.setNsearch(10);
	BCtree.setNminPnts(10);
	BCtree.setDistance(10000);
}

int iwfmMesh::findElement(double x, double y) {
	int out = -9;
	std::vector<int> id;
	BCtree.searchTree(static_cast<float>(x), static_cast<float>(y), id);

	point_b p(static_cast<float>(x), static_cast<float>(y));

	for (unsigned int i = 0; i < id.size(); ++i) {
		poly_b poly;
		int n = 4;
		if (MSH.get(id[i], 3) == 0)
			n = 3;

		for (int j = 0; j < n; ++j) {
			bg::append(poly.outer(), point_b(static_cast<float>(XY.get(MSH.get(id[i], j) - 1, 0)), static_cast<float>(XY.get(MSH.get(id[i], j) - 1, 1))));
		}
		bg::correct(poly);
		if (bg::within(p, poly)) {
			out = id[i];
			break;
		}
	}
	return out;
}

bool iwfmMesh::getPntStrat(double x, double y, std::vector<double> &elev) {
	elev.clear();
	int id = findElement(x, y);
	if (id < 0) {
		return false;
	}
	else {
		std::vector<double> wght;
		int n = 4;
		if (MSH.get(id, 3) == 0)
			n = 3;
		double sumW = 0;
		for (int i = 0; i < n; ++i) {
			double d = std::sqrt(pow(x - XY.get(MSH.get(id, i) - 1, 0), 2) + pow(y - XY.get(MSH.get(id, i) - 1, 1), 2));

			if (d < tolerance) {
				for (int j = 0; j < STRAT.ncols(); ++j) {
					elev.push_back(STRAT.get(MSH.get(id, i) - 1, j));
				}
				return true;
			}
			else {
				sumW += 1 / d;
				wght.push_back(1 / d);
			}
		}

		for (int i = 0; i < STRAT.ncols(); ++i) {
			double sumVW = 0;
			for (int j = 0; j < n; ++j) {
				sumVW += STRAT.get(MSH.get(id, j) - 1, i)*wght[j];
			}
			elev.push_back(sumVW / sumW);
		}
		return true;
	}
}

double iwfmMesh::getXY(int i, int j) {
	return XY.get(i, j);
}

int iwfmMesh::getMSH(int i, int j) {
	return MSH.get(i, j);
}

double iwfmMesh::getZ(int i, int j) {
	return STRAT.get(i, j);
}

int iwfmMesh::Nnodes() {
	return XY.nrows();
}

int iwfmMesh::Nelem() {
	return MSH.nrows();
}

int iwfmMesh::Nlay() {
	return STRAT.ncols();
}