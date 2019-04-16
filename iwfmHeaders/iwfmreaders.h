#ifndef IWFMREADERS_H
#define IWFMREADERS_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>

#include <H5Cpp.h>


struct iwfmNode{
    double X;
    double Y;
    std::vector<double> Z;
};


bool readNodeCoord(std::string filename, std::map<int, iwfmNode>& ND){
    std::ifstream datafile(filename.c_str());
    if (!datafile.good()){
        std::cout << "Can't open the file" << filename << std::endl;
        return false;
    }
    else{
        char buffer[1024];
        double xcoord, ycoord, FACT;
        int id, Nnodes;
        int cnt_lines = 0;
        while (!datafile.eof()) {
            datafile.getline(buffer,1024);
            if ('C' != buffer[0]){
                if (cnt_lines == 0){
                    std::istringstream inp(buffer);
                    inp >> Nnodes;
                    cnt_lines++;
                }
                else if (cnt_lines == 1){
                    std::istringstream inp(buffer);
                    inp >> FACT;
                    cnt_lines++;
                }
                else {
                    std::istringstream inp(buffer);
                    inp >> id;
                    inp >> xcoord;
                    inp >> ycoord;
                    iwfmNode node;
                    node.X = xcoord;
                    node.Y = ycoord;
                    ND.insert(std::pair<int, iwfmNode>(id, node));
                    cnt_lines++;
                }
            }
        }
    }
}


bool readStratigraphy(std::string filename, std::map<int, iwfmNode>& ND){
    std::ifstream datafile(filename.c_str());
    if (!datafile.good()){
        std::cout << "Can't open the file" << filename << std::endl;
        return false;
    }
    else {
        std::map<int, iwfmNode>::iterator it;
        char buffer[1024];
        double gse, a, l, FACT;
        int id, nl;
        int cnt_lines = 0;
        while (!datafile.eof()){
            datafile.getline(buffer,1024);
            if ('C' != buffer[0]){
                if (cnt_lines == 0){
                    std::istringstream inp(buffer);
                    inp >> nl;
                    cnt_lines++;
                }
                else if (cnt_lines == 1){
                    std::istringstream inp(buffer);
                    inp >> FACT;
                    cnt_lines++;
                }
                else {
                    std::istringstream inp(buffer);
                    inp >> id;
                    inp >> gse;
                    it = ND.find(id);
                    if (it != ND.end()){
                        double lay_above = gse;
                        it->second.Z.push_back(lay_above);
                        for (int i = 0; i < nl; ++i) {
                            inp >> a;
                            inp >> l;
                            lay_above = lay_above - (a + l);
                            it->second.Z.push_back(lay_above);
                        }
                    }
                }
            }
        }
    }
}

bool readMesh_Flow( std::string filename,
               std::vector<std::vector<unsigned int> >& MSH,
               std::vector<std::vector<unsigned int> >& FCELEM){
    bool out = false;

    const H5std_string FILE_NAME(filename);

    H5::H5File *GWZB_file = new H5::H5File(FILE_NAME, H5F_ACC_RDONLY);
    H5::Group *attrGroup = new H5::Group(GWZB_file->openGroup("Attributes"));
    H5::DataSet *elemDSet = new H5::DataSet(attrGroup->openDataSet("SystemData%ElementNodes"));
    H5::DataSet *faceElemDset = new H5::DataSet(attrGroup->openDataSet("SystemData%FaceElements"));

    // Read mesh and face elements
    H5::DataSpace elemSp = elemDSet->getSpace();
    H5::DataSpace facelemSp = faceElemDset->getSpace();

    // Rank is the number of dimensions of the array
    int elemRank = elemSp.getSimpleExtentNdims();
    int facelemRank = facelemSp.getSimpleExtentNdims();

    hsize_t *elem_dims = new hsize_t[elemRank];
    hsize_t *facelem_dims = new hsize_t[facelemRank];

    int ndims = elemSp.getSimpleExtentDims(elem_dims, NULL);
    ndims = facelemSp.getSimpleExtentDims(facelem_dims,NULL);
    // Allocate space for the element Node Ids
    int MSHarray[elem_dims[0]][elem_dims[1]];

    elemDSet->read(MSHarray, H5::PredType::NATIVE_INT);

    MSH.clear();
    for (unsigned int i = 0; i < static_cast<unsigned int>(elem_dims[0]); ++i) {
        std::vector<unsigned int> temp;
        for (unsigned int j = 0; j < static_cast<unsigned int>(elem_dims[1]); ++j) {
            temp.push_back(static_cast<unsigned int>(MSHarray[i][j]));
        }
        MSH.push_back(temp);
    }

    int FCELarray[facelem_dims[0]][facelem_dims[1]];
    faceElemDset->read(FCELarray,H5::PredType::NATIVE_INT);

    FCELEM.clear();
    for (unsigned int i = 0; i < static_cast<unsigned int>(facelem_dims[0]); ++i) {
        std::vector<unsigned int> temp;
        for (unsigned int j = 0; j < static_cast<unsigned int>(facelem_dims[1]); ++j) {
            temp.push_back(static_cast<unsigned int>(FCELarray[i][j]));
        }
        FCELEM.push_back(temp);
    }


    delete[] elem_dims;
    delete[] facelem_dims;
    delete  elemDSet;
    delete  faceElemDset;
    delete  attrGroup;
    delete  GWZB_file;

    out = true;
    return out;
}

bool readFaceFlows(std::string filename,
                   std::map<int, std::vector<std::vector<double> > >& FCFLOWS){
    const H5std_string FILE_NAME(filename);
    H5::H5File *GWZB_file = new H5::H5File(FILE_NAME, H5F_ACC_RDONLY);
    //Face flows [layers][time][flows between elements]
    std::vector< std::vector< std::vector< double > > > faceFlows;

    for (unsigned int i = 1; i < 5; ++i){
        std::vector< std::vector< double > > timeFlows;
        std::stringstream layName;
        layName << "Layer_" << i;
        std::cout << "Reading face flows for " << layName.str() << std::endl;
        H5::Group *laygroup = new H5::Group(GWZB_file->openGroup(layName.str()));
        H5::DataSet *layDset = new H5::DataSet(laygroup->openDataSet("FaceFlows"));

        H5::DataSpace laySp = layDset->getSpace();
        int layRank = laySp.getSimpleExtentNdims();
        hsize_t *fflow_dims = new hsize_t[layRank];
        int nlaydims = laySp.getSimpleExtentDims(fflow_dims, NULL);

        hsize_t offset[2];
        hsize_t  count[2];
        H5::DSetCreatPropList cparms = layDset->getCreatePlist();
        hsize_t chunk_dims[2];
        int rank_chunk = cparms.getChunk(2, chunk_dims);
        count[0]  = chunk_dims[0];
        count[1]  = chunk_dims[1];

        double faceFlowsArray[chunk_dims[0]][chunk_dims[1]];

        for (int it = 0; it < static_cast<int>(fflow_dims[0]); ++it){
            offset[0] = static_cast<hsize_t>(it);
            offset[1] = 0;
            laySp.selectHyperslab(H5S_SELECT_SET, count, offset);
            H5::DataSpace mspace3(rank_chunk, chunk_dims);
            layDset->read(faceFlowsArray, H5::PredType::NATIVE_DOUBLE, mspace3, laySp);
            // put all values of this time step in a vector
            std::vector<double> offsetFlows;
            for (int ix = 0; ix < static_cast<int>(chunk_dims[1]); ++ix){
                offsetFlows.push_back(faceFlowsArray[0][ix]);
            }
            timeFlows.push_back(offsetFlows);
        }
        faceFlows.push_back(timeFlows);

        delete[] fflow_dims;
        delete  layDset;
        delete laygroup;
    }

    for (int it = 0; it < faceFlows[0].size(); ++it) {
        std::vector<std::vector<double> > layflow;
        for (int ilay = 0; ilay < faceFlows.size(); ++ilay) {
            layflow.push_back(faceFlows[ilay][it]);
        }
        FCFLOWS.insert(std::pair<int, std::vector<std::vector<double> > >(it, layflow));
    }

    delete GWZB_file;

    return true;
}

void readVerticalFlows(std::string filename, std::map<int, std::vector<std::vector<double> > >& VERTFLOWS){
    const H5std_string FILE_NAME(filename);
    H5::H5File *GWZB_file = new H5::H5File(FILE_NAME, H5F_ACC_RDONLY);
    std::vector< std::vector< std::vector< double > > > vertFlows;

    for (unsigned int i = 1; i < 4; ++i){
        std::vector< std::vector< double > > timeFlows;
        std::stringstream layName;
        layName << "Layer_" << i;
        std::cout << "Reading vertical flows for " << layName.str() << std::endl;
        H5::Group *laygroup = new H5::Group(GWZB_file->openGroup(layName.str()));
        H5::DataSet *layDset = new H5::DataSet(laygroup->openDataSet("VerticalFlows"));

        H5::DataSpace laySp = layDset->getSpace();
        int layRank = laySp.getSimpleExtentNdims();
        hsize_t *vflow_dims = new hsize_t[layRank];
        int nlaydims = laySp.getSimpleExtentDims(vflow_dims, NULL);

        hsize_t offset[2];
        hsize_t  count[2];
        H5::DSetCreatPropList cparms = layDset->getCreatePlist();
        hsize_t chunk_dims[2];
        int rank_chunk = cparms.getChunk(2, chunk_dims);
        count[0]  = chunk_dims[0];
        count[1]  = chunk_dims[1];

        double vertFlowsArray[chunk_dims[0]][chunk_dims[1]];

        for (int it = 0; it < static_cast<int>(vflow_dims[0]); ++it){
            offset[0] = static_cast<hsize_t>(it);
            offset[1] = 0;
            laySp.selectHyperslab(H5S_SELECT_SET, count, offset);
            H5::DataSpace mspace3(rank_chunk, chunk_dims);
            layDset->read(vertFlowsArray, H5::PredType::NATIVE_DOUBLE, mspace3, laySp);
            // put all values of this time step in a vector
            std::vector<double> offsetFlows;
            for (int ix = 0; ix < static_cast<int>(chunk_dims[1]); ++ix){
                offsetFlows.push_back(vertFlowsArray[0][ix]);
            }
            timeFlows.push_back(offsetFlows);
        }
        vertFlows.push_back(timeFlows);

        delete[] vflow_dims;
        delete  layDset;
        delete laygroup;
    }

    for (int it = 0; it < vertFlows[0].size(); ++it) {
        std::vector<std::vector<double> > layflow;
        for (int ilay = 0; ilay < vertFlows.size(); ++ilay) {
            layflow.push_back(vertFlows[ilay][it]);
        }
        VERTFLOWS.insert(std::pair<int, std::vector<std::vector<double> > >(it, layflow));
    }

    delete GWZB_file;
}

void readDeepPercolation(std::string filename, std::map<int, std::vector<double> >& DPERC){
    const H5std_string FILE_NAME(filename);
    H5::H5File *GWZB_file = new H5::H5File(FILE_NAME, H5F_ACC_RDONLY);
    H5::Group *dpercgroup = new H5::Group(GWZB_file->openGroup("Layer_1"));
    H5::DataSet *dpercDset = new H5::DataSet(dpercgroup->openDataSet("Deep Percolation_Inflow (+)"));

    H5::DataSpace dpercSp = dpercDset->getSpace();
    int dpercRank = dpercSp.getSimpleExtentNdims();
    hsize_t *dperc_dims = new hsize_t[dpercRank];
    int ndpercdims = dpercSp.getSimpleExtentDims(dperc_dims, NULL);
    hsize_t offset[2];
    hsize_t  count[2];
    H5::DSetCreatPropList cparms = dpercDset->getCreatePlist();
    hsize_t chunk_dims[2];
    int rank_chunk = cparms.getChunk(2, chunk_dims);
    count[0]  = chunk_dims[0];
    count[1]  = chunk_dims[1];

    double dpercArray[chunk_dims[0]][chunk_dims[1]];

    for (int it = 0; it < static_cast<int>(dperc_dims[0]); ++it){
        offset[0] = static_cast<hsize_t>(it);
        offset[1] = 0;
        dpercSp.selectHyperslab(H5S_SELECT_SET, count, offset);
        H5::DataSpace mspace3(rank_chunk, chunk_dims);
        dpercDset->read(dpercArray, H5::PredType::NATIVE_DOUBLE, mspace3, dpercSp);
        std::vector<double> offsetDperc;
        for (int ix = 0; ix < static_cast<int>(chunk_dims[1]); ++ix){
            offsetDperc.push_back(dpercArray[0][ix]);
        }
        DPERC.insert(std::pair<int, std::vector<double> >(it, offsetDperc));
    }


    delete[] dperc_dims;
    delete dpercDset;
    delete dpercgroup;
    delete GWZB_file;
}

void convertNode2FaceFlow(std::map<int, std::vector<std::vector<double> > >& VERTFLOWS,
                          std::map<int, std::vector<std::vector<double> > >& VERTnodeFLOWS,
                          std::map<int, std::vector<double> >& DPERC,
                          std::vector<std::vector<unsigned int> >& MSH){

    std::map<int, std::vector<std::vector<double> > >::iterator it;
    for (it = VERTnodeFLOWS.begin(); it != VERTnodeFLOWS.end(); ++it){


    }


    std::vector<std::vector<std::vector<double> > > temp;

    for (unsigned int i = 0; i < MSH.size(); ++i) {
        std::vector<double> temp1;
        for (unsigned int j = 0; j < MSH[i].size(); ++j) {
            if (MSH[i][j] == 0)
                continue;

            int ii = MSH[i][j]-1;
            //for (int ilay = 0; ilay < VERTnodeFLOWS[ii].size(); ++ilay) {
            //    VERTnodeFLOWS[ii] =

            //}
            //VERTnodeFLOWS[ii]
            //temp.push_back()
        }
    }

}


#endif // IWFMREADERS_H
