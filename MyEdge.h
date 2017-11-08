
#ifndef _TAMAEDGE_
#define _TAMAEDGE_

#include "MyBinarySearch.h"
#include "MyVertex.h"


class MyEdge
{
  public:
    std::vector<MyVertex> v;
    std::vector<int> gvidx;
    int matid;
    
  public:
    MyEdge()
    {
    }
    MyEdge(std::vector<int> &_vidx, std::vector<MQCoordinate> &_coord, std::vector<DWORD> &_color, int _matid = -1)
    {
      Set(_vidx, _coord, _color, _matid);
    }
    
    bool operator<(const MyEdge &rhs) const
    {
      if(gvidx[0] < rhs.gvidx[0])return true;
      if(gvidx[0] == rhs.gvidx[0])
      {
        if(gvidx[1] < rhs.gvidx[1])return true;
        if(gvidx[1] == rhs.gvidx[1] && v.size() > rhs.v.size())return true;
      }
      return false;
    }
    
    void Set(std::vector<int> &_vidx, std::vector<MQCoordinate> &_coord, std::vector<DWORD> &_color, int _matid = -1)
    {
      matid = _matid;
      int N = _vidx.size();
      if(N!=_coord.size() || N!=_color.size())return;
      v.resize(N);
      for(int i=0;i<N;i++)
      {
        v[i].Set(_vidx[i], _coord[i], _color[i]);
      }
    }
    
    void Set2(std::vector<int> &_vidx, std::vector<MQCoordinate> &_coord, std::vector<DWORD> &_color, int _matid, int start, int num)
    {
      matid = _matid;
      int N = _vidx.size();
      if(N!=_coord.size() || N!=_color.size() || N<num)return;
      v.resize(num);
      int cnt = 0;
      for(int i=start;i<N && cnt<num;i++)
      {
        v[cnt].Set(_vidx[i], _coord[i], _color[i]);
        cnt++;
      }
      for(int i=0;i<N && cnt<num;i++)
      {
        v[cnt].Set(_vidx[i], _coord[i], _color[i]);
        cnt++;
      }
    }
    
    void Rotate(int newidx)
    {
      std::rotate(v.begin(),v.begin()+newidx,v.end());
    }
    
    void UpdateGVIndex(std::vector<int> &table)
    {
      int numV = v.size();
      if(numV<3)return;
      gvidx.resize(numV-1);
      int outIdx = 0;
      for(int i=0;i<numV;i++)
      {
        if(i==1)continue;
        gvidx[outIdx] = MyBinarySearch(v[i].idx, table);
        outIdx++;
      }
    }
};





#endif
