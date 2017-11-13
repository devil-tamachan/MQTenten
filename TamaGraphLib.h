
#ifndef _TAMAGRAPHLIB_
#define _TAMAGRAPHLIB_

#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include "MyGraph.h"


void SortUniqueVector(std::vector<int> &arr)
{
  std::sort(arr.begin(), arr.end());
  arr.erase(std::unique(arr.begin(), arr.end()), arr.end());
}


void MakeUniqueSortedTable_matid(std::vector<MyEdge> &edges, std::vector<int> &ret)
{
  int numEdge = edges.size();
  for(int i=0;i<numEdge;i++)
  {
    MyEdge &e = edges[i];
    ret.push_back(e.matid);
  }
  SortUniqueVector(ret);
}



void SplitGraphByMatid(std::vector<MyEdge> &edges, std::vector<std::vector<MyEdge> > &splitEdges)
{
  std::vector<int> mattable;
  MakeUniqueSortedTable_matid(edges, mattable);
  
  int numUniqueMatid = mattable.size();
  splitEdges.resize(numUniqueMatid);
  
  int numEdge = edges.size();
  for(int i=0;i<numEdge;i++)
  {
    MyEdge &e = edges[i];
    int matidx = MyBinarySearch(e.matid, mattable);
    //ASSERT(matidx!=-1);
    splitEdges[matidx].push_back(e);
  }
}
void SplitGraphByMatid(std::vector<std::vector<MyEdge> > &graphs_src, std::vector<std::vector<MyEdge> > &splitEdges)
{
  int numInGraph = graphs_src.size();
  for(int i=0;i<numInGraph;i++)
  {
    std::vector<std::vector<MyEdge> > tmpGraphs;
    SplitGraphByMatid(graphs_src[i], tmpGraphs);
    splitEdges.insert(splitEdges.end(), tmpGraphs.begin(), tmpGraphs.end());
  }
}




#endif
