#ifndef HELP_FUNC_H
#define HELP_FUNC_H

#include <boost/graph/adjacency_list.hpp>

typedef boost::adjacency_list<
        boost::vecS,
        boost::vecS,
        boost::undirectedS> graph;

void createNodeGraph(std::vector<std::vector<unsigned int> > MSH, graph& G){
    for (unsigned int i = 0; i < MSH.size(); ++i) {
        for (int j = 0; j < MSH[i].size(); ++j) {
            if (MSH[i][j] == 0)
                break;

            int jj = j+1;
            if (MSH[i][jj] == 0)
                jj = 0;
            if (j == MSH[i].size() - 1)
                jj = 0;
            boost::add_edge(MSH[i][j], MSH[i][jj], G);

        }
    }
}

std::map<int,int> list_connectedNodes(graph& G, int id){
    std::map<int,int> out;
    graph::adjacency_iterator itv, endv;
    boost::tie(itv, endv) = boost::adjacent_vertices(id, G);
    for (;itv != endv; ++ itv) {
        int icon = static_cast<int>(*itv);
        out.insert(std::pair<int,int>(icon,icon));
        std::cout << *itv << std::endl;
    }
    return out;
}

std::map<int, std::map<int,int> > ElemPerNode(std::vector<std::vector<unsigned int> >& MSH){

    std::map<int, std::map<int,int> > NDelem;
    std::map<int, std::map<int,int> >::iterator it1;


    for (int i = 0; i < MSH.size(); ++i) {
        for (int j = 0; j < MSH[i].size(); ++j) {
            if (MSH[i][j] == 0)
                break;
            it1 = NDelem.find(MSH[i][j]);
            if (it1 != NDelem.end()){
                it1->second.insert(std::pair<int,int>(i,i));
            }
            else {
                std::map<int,int> temp;
                temp.insert(std::pair<int,int>(i,i));
                NDelem.insert(std::pair<int, std::map<int,int> >(MSH[i][j], temp));
            }
        }
    }
    return NDelem;
}


#endif // HELP_FUNC_H
