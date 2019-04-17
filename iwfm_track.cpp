#include "iwfmHeaders/iwfm_track_options.h"
#include "iwfmHeaders/iwfmreaders.h"
#include "iwfmHeaders/cgal_functions.h"
#include "iwfmHeaders/tracking.h"
#include "iwfmHeaders/help_func.h"

// THis will move to a header file later
#include <osg/Geometry>
#include <osg/Geode>
#include <osgUtil/SmoothingVisitor>
#include <osgViewer/Viewer>

#include <CGAL/config.h>

int main(int argc, char *argv[]){
    iwfm_track_options ito;
    bool tf = readInputParameters(argc, argv, ito);
    if (!tf)
        return 1;


    // Read the Input files
    std::vector<std::vector<unsigned int> > MSH;
    std::vector<std::vector<unsigned int> > FCELEM;
    std::map<int, iwfmNode> ND;
    std::map<int, std::vector<std::vector<double> > > FCFLOWS;
    std::map<int, std::vector<std::vector<double> > > VERTFLOWS;
    std::map<int, std::vector<std::vector<double> > > VERTnodeFLOWS;
    std::map<int, std::vector<double> > DPERC;
    readNodeCoord(ito.Node_file, ND);
    readStratigraphy(ito.Strat_file, ND);
    readMesh_Flow(ito.GW_ZB_file, MSH, FCELEM);

    readVerticalFlows(ito.GW_ZB_file, VERTFLOWS,MSH);
    readDeepPercolation(ito.GW_ZB_file, DPERC);

    //convertNode2FaceFlow(VERTFLOWS, VERTnodeFLOWS, DPERC, MSH);

    readFaceFlows(ito.GW_ZB_file, FCFLOWS);



    std::cout << "My CGAL library is " <<  CGAL_VERSION_NR << " (1MMmmb1000)" << std::endl;
    std::cout << std::endl;
    std::cout << "where MM is the major number release, mm is the minor number release, and "
           << "b is the bug fixing number release." << std::endl;

    meshSearch mshSrch(MSH, ND, 3000);

    iwfm_PTrack iwfmtrack(mshSrch,FCELEM);
    iwfmtrack.set_particles(0,658889.816361, 4193704.600484, 50);








    return 0;




    // Visualize top and bottom layer
    osg::ref_ptr<osg::Vec3Array> vert_top = new osg::Vec3Array;//(static_cast<unsigned int>(ND.size()));
    osg::ref_ptr<osg::Vec3Array> vert_bot = new osg::Vec3Array;//(static_cast<unsigned int>(ND.size()));
    for (std::map<int, iwfmNode>::iterator it=ND.begin(); it !=ND.end(); ++it){
        vert_top->push_back(osg::Vec3(static_cast<float>(it->second.X),
                                      static_cast<float>(it->second.Y),
                                      static_cast<float>(it->second.Z[0])));
        vert_bot->push_back(osg::Vec3(static_cast<float>(it->second.X),
                                      static_cast<float>(it->second.Y),
                                      static_cast<float>(it->second.Z[4])));

        //(*vert_top)[static_cast<unsigned long>(it->first-1)].set(static_cast<float>(it->second.X),
        //                                                       static_cast<float>(it->second.Y),
        //                                                       static_cast<float>(it->second.Z[0]));
        //(*vert_bot)[static_cast<unsigned long>(it->first-1)].set(static_cast<float>(it->second.X),
        //                                                       static_cast<float>(it->second.Y),
        //                                                       static_cast<float>(it->second.Z[4]));
    }

    osg::ref_ptr<osg::DrawElementsUInt> indices = new osg::DrawElementsUInt(GL_TRIANGLES);
    for (unsigned int i = 0; i < MSH.size(); ++i) {
        if (MSH[i][3] == 0){
            (*indices).push_back(MSH[i][0]-1);
            (*indices).push_back(MSH[i][1]-1);
            (*indices).push_back(MSH[i][2]-1);
        }
        else {
            (*indices).push_back(MSH[i][0]-1);
            (*indices).push_back(MSH[i][1]-1);
            (*indices).push_back(MSH[i][2]-1);
            (*indices).push_back(MSH[i][2]-1);
            (*indices).push_back(MSH[i][3]-1);
            (*indices).push_back(MSH[i][0]-1);
        }
    }



    osg::ref_ptr<osg::Geometry> top_mesh = new osg::Geometry;
    osg::ref_ptr<osg::Geometry> bot_mesh = new osg::Geometry;
    top_mesh->setVertexArray(vert_top.get());
    top_mesh->addPrimitiveSet(indices.get());
    bot_mesh->setVertexArray(vert_bot.get());
    bot_mesh->addPrimitiveSet(indices.get());

    osgUtil::SmoothingVisitor::smooth(*top_mesh);

    osg::ref_ptr<osg::Geode> root = new osg::Geode;
    root->addDrawable(top_mesh.get());
    root->addDrawable(bot_mesh.get());

    osgViewer::Viewer viewer;
    viewer.setSceneData(root.get());
    return viewer.run();



    return 0;
}
