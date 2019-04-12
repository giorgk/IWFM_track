#ifndef CGAL_FUNCTIONS_H
#define CGAL_FUNCTIONS_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Point_set_2.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>

#include <CGAL/Barycentric_coordinates_2/Mean_value_2.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>

#include <list>

#include <Eigen/Dense>

#include "iwfmHeaders/iwfmreaders.h"



typedef CGAL::Exact_predicates_inexact_constructions_kernel                 ine_Kernel;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, ine_Kernel>   Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                            Tds;
typedef CGAL::Point_set_2<ine_Kernel,Tds>::Vertex_handle                    Vertex_handle;
typedef CGAL::Point_set_2<ine_Kernel,Tds>                                   PointSet2;
typedef ine_Kernel::Point_2                                                 ine_Point2;
typedef CGAL::Polygon_2<ine_Kernel> Polygon_2;

typedef CGAL::Barycentric_coordinates::Mean_value_2<ine_Kernel> Mean_value;
typedef CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Mean_value, ine_Kernel> Mean_value_coordinates;


void barycentricCoords(std::vector<double> xv, std::vector<double> yv, double x, double y, std::vector<double>& bcoords){
    std::vector<ine_Point2> vertices;
    for (unsigned int i = 0; i < xv.size(); ++i)
        vertices.push_back(ine_Point2(xv[i], yv[i]));


    // Instantiate the class with Wachspress coordinates for the polygon defined above.
    //Wachspress_coordinates wachspress_coordinates(vertices.begin(), vertices.end());
    //Discrete_harmonic_coordinates discrete_harmonic_coordinates(vertices.begin(), vertices.end());
    Mean_value_coordinates mean_value_coordinates(vertices.begin(), vertices.end());
    ine_Point2 p_q(x, y);
    std::vector<ine_Kernel::FT> coordinates;
    coordinates.reserve(1);

    mean_value_coordinates(p_q, std::back_inserter(coordinates),
                           CGAL::Barycentric_coordinates::UNSPECIFIED_LOCATION,
                           CGAL::Barycentric_coordinates::PRECISE);

    bcoords.clear();
    for(unsigned int j = 0; j < xv.size(); ++j){
        bcoords.push_back(CGAL::to_double(coordinates[j]));
    }
}



class meshSearch{
public:
    meshSearch(
            std::vector<std::vector<unsigned int> >& MSH_in,
            std::map<int, iwfmNode>& ND_in,
            double diam);
    int findContainerPoly(double x, double y);
    void calcUnitCoords(double x, double y, double& xu, double& yu, unsigned int elid);
    void calcBaryCoords(double x, double y, std::vector<double>& bcoords, unsigned int elid);

private:
    std::vector<std::vector<unsigned int> >& MSH;
    std::map<int, iwfmNode>& ND;
    std::vector<std::vector<double> >invCoef;

    //! A list of the barycenters of the elements
    PointSet2 BCxy;
    double search_rad;

    void prepare();
    void calculateInvCoef(std::vector<double> x, std::vector<double> y, std::vector<double>& CC);
};

meshSearch::meshSearch(std::vector<std::vector<unsigned int> >& MSH_in,
                       std::map<int, iwfmNode>& ND_in,
                       double diam)
    :
      MSH(MSH_in),
      ND(ND_in),
      search_rad(diam)
{
    prepare();
}

void meshSearch::prepare(){

    std::map<int, iwfmNode>::iterator it;
    std::vector< std::pair<ine_Point2,int> > tempxy;

    for (unsigned int i = 0; i < MSH.size(); ++i){
        double cx = 0;
        double cy = 0;
        int nc = 0;
        std::vector<double> xv, yv;
        for (unsigned int j = 0; j < MSH[i].size(); ++j){
            if (MSH[i][j] == 0)
                break;
            it = ND.find(static_cast<int>(MSH[i][j]));
            if (it != ND.end()){
                cx += it->second.X;
                cy += it->second.Y;
                nc++;
                xv.push_back(it->second.X);
                yv.push_back(it->second.Y);
            }
            else {
                std::cerr << "Node " << MSH[i][j] << " not found in the Node list" << std::endl;
            }
        }
        cx = cx/static_cast<double>(nc);
        cy = cy/static_cast<double>(nc);
        std::vector<double> CC;
        calculateInvCoef(xv, yv, CC);
        invCoef.push_back(CC);

        tempxy.push_back(std::make_pair(ine_Point2(cx, cy), i) );
    }
    BCxy.insert(tempxy.begin(), tempxy.end());

}

int meshSearch::findContainerPoly(double x, double y){

    std::list<Vertex_handle> LV;
    std::map<int, iwfmNode>::iterator ndit;

    ine_Point2 p(x, y);
    ine_Point2 p1(x - search_rad, y + search_rad);
    ine_Point2 p2(x - search_rad, y - search_rad);
    ine_Point2 p3(x + search_rad, y - search_rad);
    ine_Point2 p4(x + search_rad, y + search_rad);

    BCxy.range_search(p1, p2, p3, p4, std::back_inserter(LV));

    if (LV.size() == 0){
        return -9;
    }
    else {
        std::list<Vertex_handle>::const_iterator it = LV.begin();
        std::list<ine_Point2> pnts;
        for (;it != LV.end(); ++it){
            pnts.clear();
            int id = static_cast<int>((*it)->info());

            for (unsigned int i = 0; i < MSH[id].size(); ++i) {
                if (MSH[id][i] == 0)
                    break;
                ndit = ND.find(static_cast<int>(MSH[id][i]));
                if (ndit != ND.end()){
                    pnts.push_back(ine_Point2(ndit->second.X, ndit->second.Y));
                }
                else {
                    std::cerr << "Node " << MSH[id][i] << " not found in the Node list" << std::endl;
                }
            }

            //Polygon_2 poly(pnts.begin(),pnts.end());
            switch (CGAL::bounded_side_2(pnts.begin(), pnts.end(), p, ine_Kernel())) {
            case CGAL::ON_BOUNDED_SIDE :
                for (std::list<ine_Point2>::iterator lit = pnts.begin(); lit != pnts.end(); ++lit)
                    std::cout << lit->x() << ", " << lit->y() << std::endl;
                return id;
            case CGAL::ON_BOUNDARY:
                return id;
            case CGAL::ON_UNBOUNDED_SIDE:
                continue;
            }
        }
        return -9;
    }
}

void meshSearch::calculateInvCoef(std::vector<double> x, std::vector<double> y, std::vector<double>& CC){
    //http://www.corrmap.com/features/homography_transformation.php
    std::vector<double> Xu{0, 0, 1, 1};
    std::vector<double> Yu{0, 1, 1, 0};

    Eigen::MatrixXd A(8,8);
    Eigen::VectorXd B(8);
    Eigen::VectorXd X(8);
    for (unsigned int i = 0; i < 4; ++i){
        A(i,0) = x[i];      A(i,3) = 0;
        A(i,1) = y[i];      A(i,4) = 0;
        A(i,2) = 1;         A(i,5) = 0;
        A(i+4,3) = x[i];    A(i+4,0) = 0;
        A(i+4,4) = y[i];    A(i+4,1) = 0;
        A(i+4,5) = 1;       A(i+4,2) = 0;
        A(i,6) = -x[i]*Xu[i]; A(i,7) = -y[i]*Xu[i];
        A(i+4,6) = -x[i]*Yu[i]; A(i+4,7) = -y[i]*Yu[i];
        B(i) = Xu[i];
        B(i+4) = Yu[i];
    }

    X = A.colPivHouseholderQr().solve(B);
    //Eigen::LLT<Eigen::MatrixXd> llt;
    //llt.compute(A);
    //X = llt.solve(B);

    CC.clear();
    for (unsigned int i = 0; i < 8; ++i){
        CC.push_back(X(i));
    }
}

void meshSearch::calcUnitCoords(double x, double y, double& xu, double& yu, unsigned int elid){
    double d = invCoef[elid][6]*x + invCoef[elid][7]*y + 1;
    xu = (invCoef[elid][0]*x + invCoef[elid][1]*y +invCoef[elid][2])/d;
    yu = (invCoef[elid][3]*x + invCoef[elid][4]*y +invCoef[elid][5])/d;
}

void meshSearch::calcBaryCoords(double x, double y, std::vector<double>& bcoords, unsigned int elid){
    std::vector<double> xv;
    std::vector<double> yv;
    std::map<int, iwfmNode>::iterator it;
    int np = 4;
    if (MSH[elid][3] == 0)
        np = 3;
    for (unsigned int i = 0; i < np; ++i){
        it = ND.find(MSH[elid][i]);
        if (it != ND.end()){
            xv.push_back(it->second.X);
            yv.push_back(it->second.Y);
        }
    }
    barycentricCoords(xv, yv, x, y, bcoords);
}



#endif // CGAL_FUNCTIONS_H
