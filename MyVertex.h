
#ifndef _TAMAVERTEX_
#define _TAMAVERTEX_

class MyVertex
{
  public:
    int idx;
    MQCoordinate coord;
    DWORD color;
    
  public:
    MyVertex()
    {
    }
    MyVertex(int _idx, MQCoordinate _coord, DWORD _color)
    {
      idx = _idx;
      coord = _coord;
      color = _color;
    }
    void Set(int _idx, MQCoordinate _coord, DWORD _color)
    {
      idx = _idx;
      coord = _coord;
      color = _color;
    }
};

#endif