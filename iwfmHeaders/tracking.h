#ifndef TRACKING_H
#define TRACKING_H

//#include <Eigen/

#include "iwfmreaders.h"
#include "cgal_functions.h"

struct particle{
    double x;
    double y;
    double z;
    double v;
    double age;
    double xu;
    double yu;
    int elid;
    std::vector<double> bcxy;

};


class iwfm_PTrack{
public:
    iwfm_PTrack(meshSearch& mSrch_in,
                std::vector<std::vector<unsigned int> >& FCELEM_in);

    void track_particles();
    void set_particles( int pid,
                        double x,
                        double y,
                        double depth);
    void set_particles(std::vector<int> pid,
                       std::vector<double> x,
                       std::vector<double> y,
                       std::vector<double> Depth);

private:
    std::map<int, std::vector<particle> > streamlines;
    meshSearch& mshSrch;
    std::vector<std::vector<unsigned int> >& FCELEM;

    double convertDepth2Elev(double depth, int elId);

    void calculateBarycentricCoords(double x, double y, int el, double& xu, double& yu);

};

iwfm_PTrack::iwfm_PTrack(meshSearch& mSrch_in,
                         std::vector<std::vector<unsigned int> >& FCELEM_in)
    :
      mshSrch(mSrch_in),
      FCELEM(FCELEM)
{}



void iwfm_PTrack::set_particles(std::vector<int> pid,
                                std::vector<double> x,
                                std::vector<double> y,
                                std::vector<double> Depth){
    for (unsigned int i = 0; i < x.size(); ++i) {
        set_particles(pid[i], x[i], y[i], Depth[i]);
    }
}


void iwfm_PTrack::set_particles(int pid, double x, double y, double depth){
    int elId = mshSrch.findContainerPoly(x,y);
    if (elId>=0){
        particle p;
        p.x = x;
        p.y = y;
        double xu, yu;
        mshSrch.calcUnitCoords(x,y,xu,yu, static_cast<unsigned int>(elId));
        std::vector<double> bcoords;
        mshSrch.calcBaryCoords(x,y,bcoords, static_cast<unsigned int>(elId));
    }

    //p.z = depth;
    //p.v = 0;
    //p.age = 0;
    //p.elid = elId;

    //std::vector<particle> pv;
    //pv.push_back(p);

    //streamlines.insert(std::pair<int, std::vector<particle> >(pid, pv));
}


double iwfm_PTrack::convertDepth2Elev(double depth, int elId){

    return 0.0;
}

void iwfm_PTrack::calculateBarycentricCoords(double x, double y, int el, double& xu, double& yu){

}

#endif // TRACKING_H
