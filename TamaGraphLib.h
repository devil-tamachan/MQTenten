
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


void _ConvToMyEdge(std::vector<std::vector<int>> &src_vidx, std::vector<std::vector<MQCoordinate>> &src_coord, std::vector<std::vector<DWORD>> &src_colors, std::vector<int> &src_matids, std::vector<MyEdge> &ret_edges)
{
  //int numEdge = MIN(MIN(MIN(src_vidx.size(), src_coord.size()), src_colors.size()), src_matids.size());
  int numEdge = src_vidx.size();
  for(int i=0;i<numEdge;i++)
  {
    std::vector<int> &vidx = src_vidx[i];
    std::vector<MQCoordinate> &coord = src_coord[i];
    std::vector<DWORD> &colors = src_colors[i];
    int matid = src_matids[i];
    if(vidx.size()!=3 || coord.size()!=3 || colors.size()!=3)continue;
    MyEdge e;
    for(int k=0;k<3;k++)
    {
      MyVertex &v = e.v[k];
      v.idx = vidx[k];
      v.coord = coord[k];
      v.color = colors[k];
    }
    e.matid = matid;
    ret_edges.push_back(e);
  }
}

void SortUniqueVector(std::vector<int> &arr)
{
  std::sort(arr.begin(), arr.end());
  arr.erase(std::unique(arr.begin(), arr.end()), arr.end());
}

void MakeUniqueSortedTable_vidx(std::vector<MyEdge> &edges, std::vector<int> &ret)
{
  int numEdge = edges.size();
  for(int i=0;i<numEdge;i++)
  {
    MyEdge &e = edges[i];
    ret.push_back(e.v[0].idx);
    ret.push_back(e.v[2].idx);
  }
  SortUniqueVector(ret);
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



void SplitGraphByConnection(std::vector<MyEdge> &edges, int numV, std::vector<std::vector<MyEdge>> &splitEdges)
{
using namespace boost;
{
  typedef adjacency_list <vecS, vecS, undirectedS> Graph;
  
  Graph G(numV);
  
  int numEdge = edges.size();
  for(int i=0;i<numEdge;i++)
  {
    MyEdge &e = edges[i];
    add_edge(e.gvidx[0], e.gvidx[1], G);
  }
  
  std::vector<int> component(num_vertices(G));
  int numGraph = connected_components(G, &component[0]);
  
  splitEdges.resize(numGraph);
  numEdge = edges.size();
  for(int i=0;i<numEdge;i++)
  {
    MyEdge &e = edges[i];
    int iComp = component[e.gvidx[0]];
    splitEdges[iComp].push_back(e);
  }
}
}

void SplitGraphByMatid(std::vector<MyEdge> &edges, std::vector<std::vector<MyEdge>> &splitEdges)
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
void SplitGraphByMatid(std::vector<std::vector<MyEdge>> &graphs_src, std::vector<std::vector<MyEdge>> &splitEdges)
{
  int numInGraph = graphs_src.size();
  for(int i=0;i<numInGraph;i++)
  {
    std::vector<std::vector<MyEdge>> tmpGraphs;
    SplitGraphByMatid(graphs_src[i], tmpGraphs);
    splitEdges.insert(splitEdges.end(), tmpGraphs.begin(), tmpGraphs.end());
  }
}




#endif
