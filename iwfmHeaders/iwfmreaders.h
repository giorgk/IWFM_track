#ifndef IWFMREADERS_H
#define IWFMREADERS_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>

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


#endif // IWFMREADERS_H
