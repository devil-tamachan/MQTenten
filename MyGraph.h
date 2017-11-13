
#ifndef _TAMAGRAPH_
#define _TAMAGRAPH_


// MQTenTen
// New BSD License / GPL Licenseのデュアルライセンス。詳しくはLICENSE.txtを見てください
//
// Credit:
// toposort / dag_longest_path
// NetworkX ported from python
// networkx/networkx/algorithms/dag.py(2017/10/22) - https://github.com/networkx/networkx/blob/49fb7d68040d0709f7d12f6626b16007015d96e8/networkx/algorithms/dag.py
// Released under the 3-Clause BSD license
// Copyright (C) 2004-2017 NetworkX Developers
// Aric Hagberg <hagberg@lanl.gov>
// Dan Schult <dschult@colgate.edu>
// Pieter Swart <swart@lanl.gov>
//

#include "MyBinarySearch.h"
#include "MyVertex.h"
#include "MyEdge.h"


bool matid_comp(std::map<int,int>::value_type const &l, std::map<int,int>::value_type const &r)
{
  return l.second<r.second;
}


class MyGraph
{
  public:
    std::vector<MyEdge> edges;
    std::vector<int> table_vidx1;
    bool bPrepara;
    bool bArcSplit;
    
  public:
    MyGraph() : bPrepara(false), bArcSplit(false)
    {
    }
    MyGraph(int size) : bPrepara(false), bArcSplit(false), edges(size)
    {
    }
    
    MyEdge* new_edge()
    {
      edges.push_back(MyEdge());
      return &(edges.back());
    }
    
    void SetArcSplit(bool bEnable)
    {
      bArcSplit = bEnable;
    }
    
    void SortEdges()
    {
      std::sort(edges.begin(), edges.end());
    }
    int _MyBinarySearch(int x, int gvi)
    {
      int l  = 0;
      int r = edges.size()-1;
      
      MyEdge *_edges = &(edges[0]);
      
      while(l<r)
      {
        int c = (l+r)/2;
        if(_edges[c].gvidx[gvi] < x)
          l = c+1;
        else
          r = c;
      }
      if(_edges[l].gvidx[gvi]==x)return l;
      return -1;
    }
    int _MyBinarySearch(int x, int gvi, int l, int r)
    {
      MyEdge *_edges = &(edges[0]);
      
      while(l<r)
      {
        int c = (l+r)/2;
        if(_edges[c].gvidx[gvi] < x)
          l = c+1;
        else
          r = c;
      }
      if(_edges[l].gvidx[gvi]==x)return l;
      return -1;
    }
    void GetEdges(int u/*仮想インデックス値*/, int v/*仮想インデックス値*/, std::vector<int> &ret/*edgesのインデックス値が返る*/)
    {
      int idx0 = _MyBinarySearch(u, 0);
      if(idx0==-1)return;
      
      int rsize = 1;
      int numEdges = edges.size();
      for(int i=idx0+1;i<numEdges;i++)
      {
        MyEdge &edge = edges[i];
        if(edge.gvidx[0] != u)break;
        rsize++;
      }
      
      int idx1 = _MyBinarySearch(v, 1, idx0, idx0+rsize);
      if(idx1==-1)return;
      
      for(int i=idx1;i<numEdges;i++)
      {
        MyEdge &edge = edges[i];
        if(edge.gvidx[1] != v)break;
        ret.push_back(i);
      }
    }
    
    bool IsSimple(std::vector<int> &vidxArr)
    {
      int N = GetNumV();
      std::vector<bool> bVertex(N, false);
      int arrNum = vidxArr.size();
      for(int i=0;i<arrNum;i++)
      {
        int vidx = vidxArr[i];
        if(bVertex[vidx])return false;
        bVertex[vidx] = true;
      }
      return true;
    }
    
    int GetEdgeIdxs_Longest_Cycle(std::vector<int> &path, std::vector<int> &ret/*仮想インデックス*/, std::vector<bool> &edgesUsedCheck, int w = 1, int i = 0)
    {
      const int pathEnd = path.size()-1;
      if(i>pathEnd)return w;
      
      const int u = path[i];
      const int v = path[(i==pathEnd)? 0:i+1];
      std::vector<int> edgeIdx;
      GetEdges(u, v, edgeIdx);
      int edgeIdxSize = edgeIdx.size();
      if(edgeIdxSize<=0)
      {
        ret.clear();
        return 0;
      }
      if(edgeIdxSize==1)
      {
        int newIdx = edgeIdx[0];
        if(edgesUsedCheck[newIdx])
        {
          ret.clear();
          return 0;
        }
        edgesUsedCheck[newIdx] = true;
        ret.push_back(newIdx);
        w += edges[newIdx].v.size()-2;
        if(i+1>=pathEnd+1)
        {
          if(!IsSimple(ret))
          {
            ret.clear();
            return 0;
          }
        } else /*if(i+1<pathEnd+1)*/ {
          w = GetEdgeIdxs_Longest_Cycle(path, ret, edgesUsedCheck, w, i+1);
        }
      } else {
        if(i+1<pathEnd+1)
        {
          std::vector<int> retMax;
          std::vector<bool> edgesUsedCheckMax;
          int wMax=0;
          for(int k=0;k<edgeIdxSize;k++)
          {
            int newIdx = edgeIdx[k];
            if(edgesUsedCheck[newIdx])continue;
            std::vector<bool> edgesUsedCheckChild = edgesUsedCheck;
            edgesUsedCheckChild[newIdx] = true;
            std::vector<int> retChild = ret;
            retChild.push_back(newIdx);
            int newWeight = GetEdgeIdxs_Longest_Cycle(path, retChild, edgesUsedCheckChild, w + edges[newIdx].v.size()-2, i+1);
            
            if(newWeight>wMax)
            {
              retMax.swap(retChild);
              edgesUsedCheckMax.swap(edgesUsedCheckChild);
              wMax = newWeight;
            }
          }
          ret.swap(retMax);
          edgesUsedCheck.swap(edgesUsedCheckMax);
          w = wMax;
        } else {
          int idxMax = edgeIdx[0];
          int wMax = edges[idxMax].v.size()-2;
          for(int k=1;k<edgeIdxSize;k++)
          {
            int newIdx = edgeIdx[k];
            if(edgesUsedCheck[newIdx])continue;
            int newWeight = edges[newIdx].v.size()-2;
            if(newWeight>wMax)
            {
              idxMax = newIdx;
              wMax = newWeight;
            }
          }
          ret.push_back(idxMax);
          if(!IsSimple(ret))
          {
            ret.clear();
            return 0;
          }
          w += wMax;
        }
      }
      return w;
    }
    int GetEdgeIdxs_Longest_Arc(std::vector<int> &path, std::vector<int> &ret/*仮想インデックス*/, std::vector<bool> &edgesUsedCheck, int w = 1, int i = 0)
    {
      const int pathEnd = path.size()-1;
      if(i>pathEnd-1)return w;
      
      const int u = path[i];
      const int v = path[i+1];
      std::vector<int> edgeIdx;
      GetEdges(u, v, edgeIdx);
      int edgeIdxSize = edgeIdx.size();
      if(edgeIdxSize<=0)
      {
        ret.clear();
        return 0;
      }
      if(edgeIdxSize==1)
      {
        int newIdx = edgeIdx[0];
        if(edgesUsedCheck[newIdx])
        {
          ret.clear();
          return 0;
        }
        edgesUsedCheck[newIdx] = true;
        ret.push_back(newIdx);
        w += edges[newIdx].v.size()-2;
        if(i+1>=pathEnd)
        {
          if(!IsSimple(ret))
          {
            ret.clear();
            return 0;
          }
        } else /*if(i+1<pathEnd)*/ {
          w = GetEdgeIdxs_Longest_Arc(path, ret, edgesUsedCheck, w, i+1);
        }
      } else {
        if(i+1<pathEnd)
        {
          std::vector<int> retMax;
          std::vector<bool> edgesUsedCheckMax;
          int wMax=0;
          for(int k=0;k<edgeIdxSize;k++)
          {
            int newIdx = edgeIdx[k];
            if(edgesUsedCheck[newIdx])continue;
            std::vector<bool> edgesUsedCheckChild = edgesUsedCheck;
            edgesUsedCheckChild[newIdx] = true;
            std::vector<int> retChild = ret;
            retChild.push_back(newIdx);
            int newWeight = GetEdgeIdxs_Longest_Arc(path, retChild, edgesUsedCheckChild, w + edges[newIdx].v.size()-2, i+1);
            
            if(newWeight>wMax)
            {
              retMax.swap(retChild);
              edgesUsedCheckMax.swap(edgesUsedCheckChild);
              wMax = newWeight;
            }
          }
          ret.swap(retMax);
          edgesUsedCheck.swap(edgesUsedCheckMax);
          w = wMax;
        } else {
          int idxMax = edgeIdx[0];
          int wMax = edges[idxMax].v.size()-2;
          for(int k=1;k<edgeIdxSize;k++)
          {
            int newIdx = edgeIdx[k];
            if(edgesUsedCheck[newIdx])continue;
            int newWeight = edges[newIdx].v.size()-2;
            if(newWeight>wMax)
            {
              idxMax = newIdx;
              wMax = newWeight;
            }
          }
          ret.push_back(idxMax);
          if(!IsSimple(ret))
          {
            ret.clear();
            return 0;
          }
          w += wMax;
        }
      }
      return w;
    }
    
    struct cycle_catcher
    {
      cycle_catcher(std::vector<std::vector<int> >& _cycles/*, std::vector<std::vector<int>>& _weight*/, MyGraph *_graph, int _numV)
          : cycles(_cycles)/*, weight(_weight)*/, numMax(0), graph(_graph), numV(_numV)
      { }

      /*template <typename Path, typename Graph>
      void cycle(const Path& p, const Graph& g)
      {
        typedef typename boost::property_map<Graph, boost::vertex_index_t>::const_type IndexMap;
        IndexMap indices = boost::get(boost::vertex_index, g);
        typename Path::const_iterator i, end = p.end();
        std::vector<int> v;
        for(i = p.begin(); i != end; ++i)
        {
          int vi = get(indices, *i);
          v.push_back(vi);
        }
        int numNewV = GetWeight(v);
        if(numNewV < numMax)return;
        if(numNewV > numMax)
        {
          numMax = numNewV;
          cycles.clear();
        }
        cycles.push_back(v);
      }*/
      
      template <typename Path, typename Graph>
      void cycle(const Path& p, const Graph& g)
      {
        typedef typename boost::property_map<Graph, boost::vertex_index_t>::const_type IndexMap;
        IndexMap indices = boost::get(boost::vertex_index, g);
        typename Path::const_iterator i, end = p.end();
        std::vector<int> v;
        for(i = p.begin(); i != end; ++i)
        {
          int vi = get(indices, *i);
          v.push_back(vi);
        }
        std::vector<int> cycle2;
        std::vector<bool> edgesUsedCheck(numV, false);
        int numNewV = graph->GetEdgeIdxs_Longest_Cycle(v, cycle2, edgesUsedCheck);
        if(numNewV < numMax || numNewV==0)return;
        if(numNewV > numMax)
        {
          numMax = numNewV;
          cycles.clear();
        }
        cycles.push_back(cycle2);
      }
      /*int GetWeight(std::vector<int> &path)
      {
        int w = 1;
        int lp = path.size()-1;
        for(int i=0;i<lp;i++)
        {
          const int u = path[i];
          const int v = path[i+1];
          w += weight[u][v];
        }
        return w;
      }*/
      std::vector<std::vector<int> >& cycles;
      //std::vector<std::vector<int>>& weight;
      int numMax;
      MyGraph *graph;
      int numV;
    };
    typedef boost::directed_graph<> DiGraph;
    //typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::directedS> DiGraph;

    int _AllCycle(std::vector<std::vector<int> > &cycles, std::vector<std::vector<int> > &weight)
    {
      DiGraph g;
      
      std::vector<DiGraph::vertex_descriptor> v;
      int numV = GetNumV();
      for(int i=0;i<numV;i++)
      {
        v.push_back(g.add_vertex());
      }
      
      int numEdge = edges.size();
      for(int i=0;i<numEdge;i++)
      {
        MyEdge &edge = edges[i];
        add_edge(v[edge.gvidx[0]], v[edge.gvidx[1]], g);
      }
      
      cycle_catcher vis(cycles, this, numV);

      tiernan_all_cycles(g, vis);
      return vis.numMax;
    }
    
    
    void _BestCycle(std::vector<int> &edges2)
    {
      
      std::vector<std::vector<int> > weightmap;
      MakeWeightmap(weightmap);
      
      std::vector<std::vector<int> > cycles;
      int weight = _AllCycle(cycles, weightmap);
      
      if(cycles.size()==0)return;
      
      edges2 = cycles[0];
    }
    
    void MergeFace2_Ring(MQObject mqoOut)
    {
      SortEdges();
      std::vector<int> edges2;//edgesのインデックス値。(使用例: MyEdge &edge = edges[edges2[i]];)
      _BestCycle(edges2);
      
      if(edges2.size()>0)
      {
        int new_matid = GetMatid(edges2);
        MakeFace_Cycle(mqoOut, edges2, new_matid);
        DeleteEdges(edges2);
      }
      return;
    }
    
    void _InitAdj(std::vector<std::vector<bool> > &adj, int N)
    {
      adj.resize(N);
      for(int i=0;i<N;i++)
      {
        adj[i].resize(N, false);
      }
    }
    
    void MakeSucc(std::vector<std::vector<bool> > &succ)
    {
      int N = GetNumV();
      _InitAdj(succ, N);
      
      int numEdge = edges.size();
      for(int i=0;i<numEdge;i++)
      {
        MyEdge &e = edges[i];
        int u = e.gvidx[0];
        int v = e.gvidx[1];
        succ[u][v] = true;
      }
    }
    
    void MakePred(std::vector<std::vector<bool> > &pred)
    {
      int N = GetNumV();
      _InitAdj(pred, N);
      
      int numEdge = edges.size();
      for(int i=0;i<numEdge;i++)
      {
        MyEdge &e = edges[i];
        int u = e.gvidx[0];
        int v = e.gvidx[1];
        pred[v][u] = true;
      }
    }
    void _InitWeightmap(std::vector<std::vector<int> > &weight)
    {
      int N = GetNumV();
      std::vector<int> x(N, 1);
      weight.resize(N, x);
    }
    void MakeWeightmap(std::vector<std::vector<int> > &weight)
    {
      _InitWeightmap(weight);
      
      int numEdge = edges.size();
      for(int i=0;i<numEdge;i++)
      {
        MyEdge &e = edges[i];
        int w = e.v.size()-2;
        if(w<1)continue; //error!
        const int u = e.gvidx[0];
        const int v = e.gvidx[1];
        if(weight[u][v] < w)weight[u][v] = w;
      }
    }
    
    
    // toposort
    // NetworkX ported to C++
    // networkx/networkx/algorithms/dag.py(2017/10/22) - https://github.com/networkx/networkx/blob/49fb7d68040d0709f7d12f6626b16007015d96e8/networkx/algorithms/dag.py
    // Copyright (C) 2004-2017 NetworkX Developers
    // Aric Hagberg <hagberg@lanl.gov>
    // Dan Schult <dschult@colgate.edu>
    // Pieter Swart <swart@lanl.gov>
    void toposort(std::vector<int> &result)
    {
      int N = GetNumV();
      std::vector<bool> used(N, false);
      std::vector<int> indegree_counter(N);
      int numEdge = edges.size();
      for(int i=0;i<numEdge;i++)
      {
        MyEdge &e = edges[i];
        int v = e.gvidx[1];
        used[v] = true;
        used[e.gvidx[0]] = true;
        indegree_counter[v] = indegree_counter[v]+1;
      }
      
      std::vector<int> zero;
      for(int i=0;i<N;i++)
      {
        if(used[i]==false)continue;
        if(indegree_counter[i]==0)
        {
          zero.push_back(i);
        }
      }
      
      std::vector<std::vector<bool> > succ;
      MakeSucc(succ);
      
      while(zero.size()>0)
      {
        int node = zero.back();
        zero.pop_back();
        std::vector<bool> &succChild = succ[node];
        for(int child=0;child<N;child++)
        {
          if(succChild[child]==false)continue;
          
          int tmp = indegree_counter[child]-1;
          indegree_counter[child] = tmp;
          if(tmp==0)zero.push_back(child);
        }
        result.push_back(node);
      }
    }
    
    // dag_longest_path
    // NetworkX ported to C++
    // networkx/networkx/algorithms/dag.py(2017/10/22) - https://github.com/networkx/networkx/blob/49fb7d68040d0709f7d12f6626b16007015d96e8/networkx/algorithms/dag.py
    // Copyright (C) 2004-2017 NetworkX Developers
    // Aric Hagberg <hagberg@lanl.gov>
    // Dan Schult <dschult@colgate.edu>
    // Pieter Swart <swart@lanl.gov>
    void dag_longest_path(std::vector<int> &path)
    {
      if(edges.size()==0)return;
      
      std::vector<int> toporesult;
      toposort(toporesult);
      
      int N = GetNumV();
      
      std::vector<std::vector<bool> > pred;
      MakePred(pred);
      
      std::vector<std::vector<int> > weight;
      MakeWeightmap(weight);
      
      std::vector<std::pair<int,int> > dist(N, std::make_pair(-1,-1));
      
      int lp = toporesult.size();
      for(int i=0;i<lp;i++)
      {
        int v = toporesult[i];
        std::vector<bool> &predChild = pred[v];
        int maxu[2] = {0, v};
        for(int u=0;u<N;u++)
        {
          if(predChild[u])
          {
            int us0 = dist[u].first+weight[u][v];
            if(us0 > maxu[0])
            {
              maxu[0] = us0;
              maxu[1] = u;
            }
          }
        }
        dist[v] = (maxu[0]>=0) ? std::make_pair(maxu[0], maxu[1]) : std::make_pair(0, v);
      }
      int u = -1;
      int v = -1;
      int maxcost = -1;
      for(int i=0;i<N;i++)
      {
        std::pair<int,int> &d = dist[i];
        if(d.first>maxcost)
        {
          maxcost = d.first;
          v = i;
        }
      }
      while(u!=v)
      {
        path.push_back(v);
        u=v;
        v=dist[v].second;
      }
      std::reverse(path.begin(), path.end());
    }
    
    void MergeFace2_Arc(MQObject mqoOut)
    {
      std::vector<int> edges2;
      dag_longest_path(edges2);
      
      int numV = GetNumV();
      std::vector<int> arcEdgeIdx;
      std::vector<bool> edgesUsedCheck(numV, false);
      int numNewV = GetEdgeIdxs_Longest_Arc(edges2, arcEdgeIdx, edgesUsedCheck);
      if(numNewV==0)return;
      
      if(arcEdgeIdx.size()>0)
      {
        int new_matid = GetMatid(arcEdgeIdx);
        MakeFace_Arc(mqoOut, arcEdgeIdx, new_matid);
        DeleteEdges(arcEdgeIdx);
      }
      return;
    }
    
    int GetMatid(std::vector<int> &idxEdges)
    {
      std::map<int,int> matids;
      int N = idxEdges.size();
      if(N==0)return -1;
      for(int i=0;i<N;i++)
      {
        int matid = edges[idxEdges[i]].matid;
        int v = 0;
        if(matids.find(matid) != matids.end())
          v = matids[matid];
        matids[matid] = v;
      }
      std::map<int,int>::iterator it = std::max_element(matids.begin(), matids.end(), matid_comp);
      return it->first;
    }
    
    void DeleteEdges(std::vector<int> &idxEdges)
    {
      std::sort(idxEdges.begin(), idxEdges.end());
      int numDel = idxEdges.size();
      for(int i=numDel-1;i>=0;i--)
      {
        edges.erase(edges.begin() + idxEdges[i]);
      }
    }
    
    void _MakeFace_OutVertex(MyVertex &v, std::vector<int> &vidx, std::vector<MQCoordinate> &coord, std::vector<DWORD> &colors)
    {
      vidx.push_back(v.idx);
      coord.push_back(v.coord);
      colors.push_back(v.color);
    }
    
    void _MakeFace_ReadEdges(std::vector<int> &idxEdges, std::vector<int> &vidx, std::vector<MQCoordinate> &coord, std::vector<DWORD> &colors, bool bAcyclic = false)
    {
      int idxEdgesSize = idxEdges.size();
      for(int i=0;i<idxEdgesSize;i++)
      {
        MyEdge &e = edges[idxEdges[i]];
        _MakeFace_OutVertex(e.v[0], vidx, coord, colors);
        for(int k=e.v.size()-1;k>=3;k--) _MakeFace_OutVertex(e.v[k], vidx, coord, colors);
        if(i==idxEdgesSize-1 && bAcyclic)_MakeFace_OutVertex(e.v[2], vidx, coord, colors);
      }
      std::reverse(vidx.begin(),vidx.end());
      std::reverse(coord.begin(),coord.end());
      std::reverse(colors.begin(),colors.end());
    }
    
    void _MakeFace_ToMetaseq(MQObject mqoOut, std::vector<int> &vidx, std::vector<MQCoordinate> &coord, std::vector<DWORD> &colors, int new_matid, int start, int num)
    {
      if(vidx.size()<3 || num<3)return;
      int addfi = mqoOut->AddFace(num<=0?vidx.size():num, &(vidx[start]));
      mqoOut->SetFaceCoordinateArray(addfi, &(coord[start]));
      SetVertexColors2(mqoOut, addfi, colors, start, num<=0?colors.size():num);
      if(addfi!=-1 && new_matid!=-1)mqoOut->SetFaceMaterial(addfi, new_matid);
    }
    void _MakeFace_ToMetaseq(MQObject mqoOut, std::vector<int> &vidx, std::vector<MQCoordinate> &coord, std::vector<DWORD> &colors, int new_matid)
    {
      _MakeFace_ToMetaseq(mqoOut, vidx, coord, colors, new_matid, 0, vidx.size());
    }
    
    bool _MakeFace_IsSplit(MQObject mqoOut, std::vector<MQPoint> &pts, int center, MQPoint delPts)
    {
      MQPoint &ptsLeft = pts[center-1];
      MQPoint &ptsCheck= pts[center];
      MQPoint &ptsRight= pts[center+1];
      MQPoint n1 = GetNormal(ptsLeft, ptsCheck, ptsRight);
      MQPoint n2 = GetNormal(ptsLeft, delPts, ptsRight);
      double size = GetSize(n1) * GetSize(n2);
      if(size<=0.0)return true;
      double c = GetInnerProduct(n1, n2) / size;
      if(c>0.5)return true;
      return false;
    }
    void _MakeFace_ToMetaseq_Arc_EnableSplit(MQObject mqoOut, std::vector<int> &vidx, std::vector<MQCoordinate> &coord, std::vector<DWORD> &colors, int new_matid, MQPoint delPts)
    {
      int l=0;
      int numVidx = vidx.size();
      if(numVidx<3)return;
      std::vector<MQPoint> pts;
      GetPointByVidx(mqoOut, vidx, pts);
      for(int center=1;center<numVidx-1;center++)
      {
        if(_MakeFace_IsSplit(mqoOut, pts, center, delPts))
        {
          int outSize = center-l+1;
          _MakeFace_ToMetaseq(mqoOut, vidx, coord, colors, new_matid, l, outSize);
          l = center;
        }
      }
      _MakeFace_ToMetaseq(mqoOut, vidx, coord, colors, new_matid, l, numVidx-l);
    }
    void _MakeFace_ToMetaseq_Arc(MQObject mqoOut, std::vector<int> &vidx, std::vector<MQCoordinate> &coord, std::vector<DWORD> &colors, int new_matid, MQPoint delPts)
    {
      if(bArcSplit)_MakeFace_ToMetaseq_Arc_EnableSplit(mqoOut, vidx, coord, colors, new_matid, delPts);
      else         _MakeFace_ToMetaseq                (mqoOut, vidx, coord, colors, new_matid);
    }
    
    void MakeFace_Cycle(MQObject mqoOut, std::vector<int> &idxEdges, int new_matid)
    {
      std::vector<int> vidx;
      std::vector<MQCoordinate> coord;
      std::vector<DWORD> colors;
      
      _MakeFace_ReadEdges(idxEdges, vidx, coord, colors, false);
      _MakeFace_ToMetaseq(mqoOut, vidx, coord, colors, new_matid);
    }
    
    void MakeFace_Arc(MQObject mqoOut, std::vector<int> &idxEdges, int new_matid)
    {
      std::vector<int> vidx;
      std::vector<MQCoordinate> coord;
      std::vector<DWORD> colors;
      
      if(idxEdges.size()<=0)return;
      MyEdge &e = edges[idxEdges[0]];
      int delVi = e.v[1].idx;
      MQPoint delPts = GetPointByVi(mqoOut, delVi);
      _MakeFace_ReadEdges(idxEdges, vidx, coord, colors, true);
      _MakeFace_ToMetaseq_Arc(mqoOut, vidx, coord, colors, new_matid, delPts);
    }
    
    int GetNumV()
    {
      if(!bPrepara)Prepara();
      return table_vidx1.size();
    }
    
    void Prepara()
    {
      _MakeUniqueSortedTable_vidx();
      _UpdateGVIndex();
      bPrepara = true;
    }
    void SplitByConnection(std::vector<MyGraph> &splitGraphs)
    {
      using namespace boost;
      
      typedef adjacency_list <vecS, vecS, undirectedS> Graph;
      
      
      int numEdge = edges.size();
      if(numEdge==0)return;
      
      Graph G(GetNumV());
      for(int i=0;i<numEdge;i++)
      {
        MyEdge &e = edges[i];
        add_edge(e.gvidx[0], e.gvidx[1], G);
      }
      
      std::vector<int> component(num_vertices(G));
      int numGraph = connected_components(G, &component[0]);
      
      splitGraphs.resize(numGraph);
      numEdge = edges.size();
      for(int i=0;i<numEdge;i++)
      {
        MyEdge &e = edges[i];
        int iComp = component[e.gvidx[0]];
        splitGraphs[iComp].edges.push_back(e);
      }
      for(int i=0;i<numGraph;i++)splitGraphs[i].SetArcSplit(bArcSplit);
    }
    void SplitByMatid(std::vector<MyGraph> &splitGraphs)
    {
      std::vector<int> mattable;
      _MakeUniqueSortedTable_matid(mattable);
      
      int numUniqueMatid = mattable.size();
      splitGraphs.resize(numUniqueMatid);
      
      int numEdge = edges.size();
      for(int i=0;i<numEdge;i++)
      {
        MyEdge &e = edges[i];
        int matidx = MyBinarySearch(e.matid, mattable);
        //ASSERT(matidx!=-1);
        splitGraphs[matidx].edges.push_back(e);
      }
      for(int i=0;i<numUniqueMatid;i++)splitGraphs[i].SetArcSplit(bArcSplit);
    }
    
  private:
    static void _SortUniqueVector(std::vector<int> &arr)
    {
      std::sort(arr.begin(), arr.end());
      arr.erase(std::unique(arr.begin(), arr.end()), arr.end());
    }
    void _MakeUniqueSortedTable_vidx()
    {
      int numEdge = edges.size();
      for(int i=0;i<numEdge;i++)
      {
        MyEdge &e = edges[i];
        int numV = e.v.size();
        for(int k=0;k<numV;k++)
        {
          if(k==1)continue;
          table_vidx1.push_back(e.v[k].idx);
        }
      }
      _SortUniqueVector(table_vidx1);
    }
    void _MakeUniqueSortedTable_matid(std::vector<int> &ret)
    {
      int numEdge = edges.size();
      for(int i=0;i<numEdge;i++)
      {
        MyEdge &e = edges[i];
        ret.push_back(e.matid);
      }
      _SortUniqueVector(ret);
    }
    
    void _UpdateGVIndex()
    {
      int num = edges.size();
      for(int i=0;i<num;i++)
      {
        edges[i].UpdateGVIndex(table_vidx1);
      }
    }
};








#endif
