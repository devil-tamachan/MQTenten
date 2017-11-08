

//　問題点もりだくさん
//
// 塊が２個とかあるとおかしくなる。選択面を赤く表示とかもうちょっとなんとかしないと
// そもそも面がつながってる塊検知するのどうすりゃいいの？
//
// 真ん中に穴、コの字型、頂点から４個以上のメッシュがでてるとかだと分割がうまくいかない
// ↑Constrainedで境界エッジを強制的に作るにしても穴の検知がめんどくさそう
// 
// UVの問題点も盛り沢山
// そもそも同じ頂点でUV値が違う面が複数繋がってる場合は、境界をConstrainedに追加して強制的にエッジを生成しないといけない
//
// UVを調べる下のアルゴリズムは不完全。ADEがBCEになる可能性もある。
//    面の法線、E->面中心点の角度を生成して、面の法線が一定角度以内（+-45度とか？これで反転ポリゴン除外・・・そもそも表ポリゴンしかないから必要ないか・・・）、一番近い角度のUVを当てるみたいなアルゴリズムがいいかも
// A   B  C   D
//
//       E

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include "MQBasePlugin.h"
#include "MQWidget.h"
#include "MQ3dLib.h"
#include "TamaMQLib.h"

BOOL TenTen(MQDocument doc);


#include <iostream>
#include <list>

#ifdef ENABLE_TENUCHI

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Simple_cartesian<double> K;

typedef K::FT FT;
typedef K::Line_3 Line3;
typedef K::Point_2 Point2;
typedef K::Point_3 Point3;
typedef K::Triangle_3 Triangle;
typedef K::Direction_3 Direction;
typedef K::Vector_3 Vector3;

#endif //#ifdef ENABLE_TENUCHI

#include <boost/graph/directed_graph.hpp>
#include <boost/graph/tiernan_all_cycles.hpp>

class MQPINFO
{
public:
  MQPINFO()
  {
    i = -1;
    oidx = -1;
  }
  MQPINFO(int _oidx, int _i)
  {
    i = _i;
    oidx = _oidx;
  }
  int i;
  int oidx;
};

#ifdef ENABLE_TENUCHI
typedef CGAL::Triangulation_vertex_base_with_info_2<MQPINFO, K> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> DT;
#endif //#ifdef ENABLE_TENUCHI

#include "TamaGraphLib.h"


BOOL APIENTRY DllMain(HANDLE hModule, DWORD ul_reason_for_call, LPVOID lpReserved)
{
	//プラグインとしては特に必要な処理はないので、何もせずにTRUEを返す
    return TRUE;
}

class TenTenWindow : public MQWindow
{
public:
  TenTenWindow(MQWindowBase& parent);
  ~TenTenWindow();
  
  MQButton *addvbtn;
  MQButton *delvbtn;
  MQButton *reEdgebtn;
  MQCheckBox *chkArcSplit;
  int mode;
  bool bArcSplit;
  
  BOOL OnClickAddV(MQWidgetBase *sender, MQDocument doc)
  {
    addvbtn->SetDown(true);
    delvbtn->SetDown(false);
    reEdgebtn->SetDown(false);
    mode=0;
    return FALSE;
  }
  BOOL OnClickDelV(MQWidgetBase *sender, MQDocument doc)
  {
    delvbtn->SetDown(true);
#ifdef ENABLE_TENUCHI
    addvbtn->SetDown(false);
    reEdgebtn->SetDown(false);
#endif //#ifdef ENABLE_TENUCHI
    mode=1;
    return FALSE;
  }
  BOOL OnClickReEdge(MQWidgetBase *sender, MQDocument doc)
  {
    addvbtn->SetDown(false);
    delvbtn->SetDown(false);
    reEdgebtn->SetDown(true);
    mode=2;
    return FALSE;
  }
  BOOL OnCheckArcSplit(MQWidgetBase *sender, MQDocument doc)
  {
    bArcSplit = !bArcSplit;
    chkArcSplit->SetChecked(bArcSplit);
    return FALSE;
  }
};

TenTenWindow::TenTenWindow(MQWindowBase& parent) : MQWindow(parent), addvbtn(NULL), delvbtn(NULL), reEdgebtn(NULL), mode(0), bArcSplit(false)
{
#ifdef ENABLE_TENUCHI
  mode=0;
#else
  mode=1;
#endif //#ifdef ENABLE_TENUCHI
  SetTitle(L"TenTen");

  MQFrame *mainFrame = CreateHorizontalFrame(this);

  MQFrame *paramFrame = CreateVerticalFrame(mainFrame);
  //paramFrame->SetMatrixColumn(2);

#ifdef ENABLE_TENUCHI
  addvbtn = CreateButton(paramFrame, L"点打ち");
  addvbtn->SetToggle(true);
  addvbtn->SetDown(true);
#endif //#ifdef ENABLE_TENUCHI
  delvbtn = CreateButton(paramFrame, L"点消し");
  delvbtn->SetToggle(true);
#ifndef ENABLE_TENUCHI
  delvbtn->SetDown(true);
#endif
#ifdef ENABLE_TENUCHI
  reEdgebtn = CreateButton(paramFrame, L"辺貼り直し");
  reEdgebtn->SetToggle(true);
#endif //#ifdef ENABLE_TENUCHI
  
#ifdef ENABLE_TENUCHI
  addvbtn->AddClickEvent(this, &TenTenWindow::OnClickAddV);
  reEdgebtn->AddClickEvent(this, &TenTenWindow::OnClickReEdge);
#endif //#ifdef ENABLE_TENUCHI
  delvbtn->AddClickEvent(this, &TenTenWindow::OnClickDelV);
  
  chkArcSplit = CreateCheckBox(paramFrame, L"弧分割");
  chkArcSplit->SetChecked(bArcSplit);
  chkArcSplit->AddChangedEvent(this, &TenTenWindow::OnCheckArcSplit);

}

TenTenWindow::~TenTenWindow()
{
}

class TenTenPlugin : public MQCommandPlugin
{
public:
  TenTenPlugin() : bActive(false), win(NULL), errMsg(0)
  {
  }
  
  virtual void GetPlugInID(DWORD *Product, DWORD *ID);
  virtual const char *GetPlugInName(void);
  virtual const char *EnumString(void);
  virtual BOOL Initialize();
  virtual void Exit();
  virtual BOOL Activate(MQDocument doc, BOOL flag);
  virtual void OnDraw(MQDocument doc, MQScene scene, int width, int height);
  
  BOOL AddVertexOnTri(MQDocument doc, int oi, int fi, MQPoint &hitpoint, bool bUpdateUndo);
  
  // 左ボタンが押されたとき
  virtual BOOL OnLeftButtonDown(MQDocument doc, MQScene scene, MOUSE_BUTTON_STATE& state);
  // 左ボタンが押されながらマウスが移動したとき
  //virtual BOOL OnLeftButtonMove(MQDocument doc, MQScene scene, MOUSE_BUTTON_STATE& state);
  // 左ボタンが離されたとき
  //virtual BOOL OnLeftButtonUp(MQDocument doc, MQScene scene, MOUSE_BUTTON_STATE& state);
  
  int GetMode()
  {
    if(win)return win->mode;
    else return -1;
  }
  int GetEnableArcSplit()
  {
    if(win)return win->bArcSplit;
    else return false;
  }
  
private:
  BOOL bActive;
  TenTenWindow *win;
  int errMsg;

};

void TenTenPlugin::GetPlugInID(DWORD *Product, DWORD *ID)
{
  // プロダクト名(制作者名)とIDを、全部で64bitの値として返す
  // 値は他と重複しないようなランダムなもので良い
  *Product = 0xA8BEE201;
  *ID      = 0x9A9D0493;
}

const char *TenTenPlugin::GetPlugInName(void)
{
  // プラグイン名
  return "MQTenTen           Copyright(C) 2017, tamachan";
}

const char *TenTenPlugin::EnumString()
{
  return "TenTen";
}

BOOL TenTenPlugin::Initialize()
{
  // 特に何もしないので、そのままTRUEを返す
  return TRUE;
}

void TenTenPlugin::Exit()
{
  if(win)
  {
    delete win;
    win = NULL;
  }
}

MQBasePlugin *GetPluginClass()
{
	static TenTenPlugin plugin;
	return &plugin;
}
/*
void AABB_AddMQObj2(std::list<Triangle> &triangles, MQDocument doc, MQObject o, int oi)
{
  if(o==NULL)return;
  int numV = o->GetVertexCount();
  int numF = o->GetFaceCount();
  for(int k=0;k<numF;k++)
  {
    if(doc->IsSelectFace(oi, k)==FALSE)continue;
    int numFV = o->GetFacePointCount(k);
    if(numFV<3)continue;
    
    std::vector<int> index(numFV);
    o->GetFacePointArray(k, &(*index.begin()));
    
    int numTri = numFV - 2;
    std::vector<int> indices(numTri*3);
    if(numFV==3)
    {
      indices[0] = 0;
      indices[1] = 1;
      indices[2] = 2;
    } else {
      std::vector<MQPoint> pts(numFV);
      for(int i=0; i<numFV; i++)
      {
        pts[i] = o->GetVertex(index[i]);
      }
      doc->Triangulate(&(*pts.begin()), numFV, &(*indices.begin()), numTri*3);
    }
    for(int m=0;m<numTri;m++)
    {
      MQPoint p;
      p = o->GetVertex(index[indices[m*3+0]]);
      Point a(p.x, p.y, p.z);
      p = o->GetVertex(index[indices[m*3+1]]);
      Point b(p.x, p.y, p.z);
      p = o->GetVertex(index[indices[m*3+2]]);
      Point c(p.x, p.y, p.z);
      
      triangles.push_back(Triangle(a,b,c));
    }
  }
}

void AABB_AddMQObjs(std::list<Triangle> &triangles, MQDocument doc)
{
  //std::list<Triangle> triangles;
  int numobj = doc->GetObjectCount();
  for(int oi=0;oi<numobj;oi++)
  {
    MQObject o = doc->GetObject(oi);
    if(o==NULL || o->GetLocking() || o->GetVisible()==0)continue;
    AABB_AddMQObj2(triangles, doc, o, oi);
  }
  Tree tree(triangles.begin(),triangles.end());
}
*/

bool IsVisibleFace(MQScene scene, MQObject o, int fi)
{
  int numFV = o->GetFacePointCount(fi);
  std::vector<int> vidx(numFV);
  o->GetFacePointArray(fi, &(*vidx.begin()));

  std::vector<MQPoint> p2d(numFV);
  for(int i=0;i<numFV;i++)
  {
    p2d[i] = scene->Convert3DToScreen(o->GetVertex(vidx[i]));
    if(p2d[i].z <= 0)return false;
  }

  return ((p2d[0].x - p2d[2].x) * (p2d[1].y - p2d[2].y) - (p2d[1].x - p2d[2].x) * (p2d[0].y - p2d[2].y)) > 0;
}

bool IsFrontTriSelected(MQDocument doc, MQScene scene)
{
  int numobj = doc->GetObjectCount();
  for(int oi=0;oi<numobj;oi++)
  {
    MQObject o = doc->GetObject(oi);
    if(o==NULL || o->GetLocking() || o->GetVisible()==0)continue;
    
    int numF = o->GetFaceCount();
    for(int fi=0;fi<numF;fi++)
    {
      if(doc->IsSelectFace(oi, fi)==FALSE)continue;
      int numFV = o->GetFacePointCount(fi);
      if(numFV==3 && IsVisibleFace(scene, o, fi)==false)return false;
    }
  }
  return true;
}



//---------------------------------------------------------------------------
//  TenTenPlugin::Activate
//    表示・非表示切り替え要求
//---------------------------------------------------------------------------
BOOL TenTenPlugin::Activate(MQDocument doc, BOOL flag)
{
  if(flag)
  {
    if(win==NULL)win = new TenTenWindow(MQWindow::GetMainWindow());
    if(win)win->SetVisible(true);
  } else {
    if(win)win->SetVisible(false);
  }
  bActive = flag;
  RedrawAllScene();
  /*
  if (!flag){
    if(m_bDragging){
      m_bDragging = false;
      RedrawAllScene();
    }
  }*/
  
  // そのままflagを返す
  return flag;
}

void TenTenPlugin::OnDraw(MQDocument doc, MQScene scene, int width, int height)
{
  if (bActive && errMsg)
  {
    DRAWING_TEXT_PARAM param;
    param.HorzAlign = TEXT_ALIGN_CENTER;
    param.VertAlign = TEXT_ALIGN_CENTER;
    param.Color = MQColor(0,0,0);
    param.FontScale = 1.0f;
    param.ScreenPos = MQPoint(width/2.0, height/2.0, 0);
    CreateDrawingText(doc, L"選択面に裏面ポリゴンがあります", param);
    param.Color = MQColor(1,1,0);
    param.FontScale = 1.1f;
    CreateDrawingText(doc, L"選択面に裏面ポリゴンがあります", param);
    errMsg = 0;
  }
}

MQCoordinate CalcCoord(MQPoint &p1, MQPoint &p2, MQPoint &p3, MQPoint &pn, std::vector<MQCoordinate> &uv)
{
  MQCoordinate ret;
  
  float la = GetSize(GetCrossProduct(p1-p3, p1-p2));
  
  if(la==0.0f)
  {
    ret.u = uv[0].u;
    ret.v = uv[0].v;
    return ret;
  }
  
  MQPoint v1 = p1-pn;
  MQPoint v2 = p2-pn;
  MQPoint v3 = p3-pn;
  
  float l3 = GetSize(GetCrossProduct(v1, v2)) / la;
  float l1 = GetSize(GetCrossProduct(v2, v3)) / la;
  float l2 = GetSize(GetCrossProduct(v3, v1)) / la;
  
  ret.u = uv[0].u * l1 + uv[1].u * l2 + uv[2].u * l3;
  ret.v = uv[0].v * l1 + uv[1].v * l2 + uv[2].v * l3;
  
  return ret;
}

BOOL TenTenPlugin::AddVertexOnTri(MQDocument doc, int oi, int fi, MQPoint &hitpoint, bool bUpdateUndo = true)
{
    MQObject o = doc->GetObject(oi);
    if(o==NULL || doc->IsSelectFace(oi, fi)==0)return FALSE;
    int idx[3];
    
    int numfv = o->GetFacePointCount(fi);
    if(numfv!=3)return FALSE;
    
    idx[2] = o->AddVertex(hitpoint);
    std::vector<int> vidxOld(numfv+1);
    o->GetFacePointArray(fi, &(*vidxOld.begin()));
    std::vector<MQCoordinate> vcoordOld(numfv+1);
    o->GetFaceCoordinateArray(fi, &(*vcoordOld.begin()));
    vidxOld[numfv] = vidxOld[0];
    vcoordOld[numfv] = vcoordOld[0];
    MQCoordinate coord[3];
    coord[2] = CalcCoord(o->GetVertex(vidxOld[0]), o->GetVertex(vidxOld[1]), o->GetVertex(vidxOld[2]), hitpoint, vcoordOld);
    int matidx = o->GetFaceMaterial(fi);
    for(int i=0;i<numfv;i++)
    {
      idx[0] = vidxOld[i];
      idx[1] = vidxOld[i+1];
      int fi = o->AddFace(3, idx);
      doc->AddSelectFace(oi, fi);
      coord[0] = vcoordOld[i];
      coord[1] = vcoordOld[i+1];
      o->SetFaceCoordinateArray(fi, coord);
      o->SetFaceMaterial(fi, matidx);
    }
    for(int i=0;i<numfv;i++)doc->AddSelectVertex(oi, idx[i]);
    doc->AddSelectVertex(oi, idx[2]);
    
    o->DeleteFace(fi);

    //doc->Compact();
    
    RedrawAllScene();
    if(bUpdateUndo)
    {
      OutputDebugStringA("Undo10\n");
      UpdateUndo();
    }
    
    return TRUE;
}

class CReEdgeVertex
{
public:
  int vidx;
  MQPoint p;
  MQObject o;
  
  CReEdgeVertex() : vidx(-1)
  {
  }
  
  CReEdgeVertex(int _vidx)
  {
    vidx = _vidx;
  }
  CReEdgeVertex(int _vidx, MQObject _o, MQPoint _p)
  {
    vidx = _vidx;
   p = _p;
   o = _o;
  }
  
  bool operator < (const CReEdgeVertex& r) const
  {
    return vidx < r.vidx;
  }
  bool operator == (const CReEdgeVertex& r) const
  {
    return vidx == r.vidx;
  }
};

int RightIndex(int i, int num)
{
  int r = i+1;
  if(r>=num)r=0;
  return r;
}
int LeftIndex(int i, int num)
{
  int left = i-1;
  if(left<0)left=num-1;
  return left;
}

int _MatchIndexLR(int idxl, int idxc, int idxr, std::vector<int> &vidx, int *reti)
{
  *reti = -1;
  int num = vidx.size();
  if(num<=1)return 0;
  for(int i=0;i<num;i++)
  {
    if(vidx[i]==idxc)
    {
      int l = i-1;
      if(l<0)l = num-1;
      int r = i+1;
      if(r==num)r=0;

      int ret = 2;
      ret += (vidx[r]==idxr) ? 0x1 : 0;
      ret += (vidx[l]==idxl) ? 0x4 : 0;
      *reti = i;

      return ret;
    }
  }
  return 0;
}

class FoundUV
{
public:
  MQCoordinate uv;
  int vidx;
  bool bFound;
  int bestmatch;

  FoundUV() : bFound(false), vidx(-1), bestmatch(0)
  {
  }
};


void _FindBestUV(MQObject o, std::vector<FoundUV> &retUV, int idxc_retUV)
{
  int numRetUV = retUV.size();
  int idxc = retUV[idxc_retUV].vidx;
  
  int idxr_retUV = RightIndex(idxc_retUV, numRetUV);
  int idxr = retUV[idxr_retUV].vidx;
  
  int idxl_retUV = LeftIndex(idxc_retUV, numRetUV);
  int idxl = retUV[idxl_retUV].vidx;

  int numf = o->GetVertexRelatedFaces(idxc, NULL);
  std::vector<int> _fidx(numf);
  o->GetVertexRelatedFaces(idxc, &(*_fidx.begin()));
  std::vector<int> fidx(numf);

  //_fidxを三角形以上に選別
  for(int i=0;i<numf;i++)
  {
    if(o->GetFacePointCount(_fidx[i])>=3)
      fidx.push_back(_fidx[i]);
  }
  numf = fidx.size();

  std::vector<std::vector<int>> vidxArr(numf);
  
  int besti = 0;
  int bestmatch = 0;
  for(int i=0;i<numf;i++)
  {
    int fi = fidx[i];
    std::vector<int> &vidx = vidxArr[i];
    int numV = o->GetFacePointCount(fi);
    vidx.resize(numV);
    o->GetFacePointArray(fi, &(*vidx.begin()));

    int reti = -1;
    int match = _MatchIndexLR(idxl, idxc, idxr, vidx, &reti);
    switch(match)
    {
    case 0: //見つからない。来ないはず
    default://来ないはず
    case 1: //rのみマッチ。来ないはず
    case 4: //lのみマッチ。来ないはず
    case 5: //l+rマッチ。来ないはず
      continue;
    case 2: //cのみマッチ
    case 3: //c+rマッチ
    case 6: //l+cマッチ
      if(match>bestmatch)
      {
        bestmatch = match;
        besti = i;
      }
      continue;
    case 7: //l+c+rマッチ
      bestmatch = 7;
      besti = i;
      goto exit__FindBestUV;
    }
  }
exit__FindBestUV:
  if(bestmatch>0)
  {

  }
  //対象頂点の右と左頂点が一緒
  //対象頂点の左が一緒
  //対象頂点の右が一緒
  //適当にはじめに見つかった頂点のUV
}
/*
MQCoordinate FindBestUV(MQObject o, std::vector<FoundUV> &retUV)
{
  int numV = retUV.size();
  for(int i=0;i<numV;i++)
  {
    FoundUV &tv = retUV[i];
    if(tv.bFound==false)_FindBestUV(o, retUV, i);
  }

}*/

#ifdef ENABLE_TENUCHI
BOOL ReEdge(MQDocument doc, MQScene scene)
{
  int numobj = doc->GetObjectCount();
  std::vector<std::vector<CReEdgeVertex>> vidx(numobj);
  for(int oi=0;oi<numobj;oi++)
  {
    MQObject o = doc->GetObject(oi);
    if(o==NULL || o->GetLocking() || o->GetVisible()==0)continue;
    
    int numF = o->GetFaceCount();
    for(int fi=0;fi<numF;fi++)
    {
      if(doc->IsSelectFace(oi, fi)==FALSE)continue;
      int numFV = o->GetFacePointCount(fi);
      std::vector<int> vidx2(numFV);
      o->GetFacePointArray(fi, &(*vidx2.begin()));
      for(int k=0;k<vidx2.size();k++)
      {
        int t = vidx2[k];
        vidx[oi].push_back(CReEdgeVertex(t, o, o->GetVertex(t)));
      }
      o->DeleteFace(fi, false);
    }
  }
  
  for(int oi=0;oi<numobj;oi++)
  {
    std::sort(vidx[oi].begin(), vidx[oi].end());
    vidx[oi].erase(std::unique(vidx[oi].begin(), vidx[oi].end()), vidx[oi].end());
  }
  std::vector<std::pair<Point2,MQPINFO>> points;
  for(int oi=0;oi<numobj;oi++)
  {
    int numV = vidx[oi].size();
    for(int i = 0;i<numV;i++)
    {
      MQPoint p2d = scene->Convert3DToScreen(vidx[oi][i].p);
      points.push_back(std::make_pair(Point2(p2d.x, p2d.y),MQPINFO(oi, i)));
    }
  }
  DT dt;
  dt.insert(points.begin(), points.end());
  DT::Finite_faces_iterator fit;
  for (fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); ++fit)
  {
    int idx[3];
    DT::Face_handle face = fit;
    MQPINFO info0 = face->vertex(0)->info();
    MQPINFO info2 = face->vertex(2)->info();
    MQPINFO info1 = face->vertex(1)->info();
    CReEdgeVertex &rv0 = vidx[info0.oidx][info0.i];
    CReEdgeVertex &rv1 = vidx[info1.oidx][info1.i];
    CReEdgeVertex &rv2 = vidx[info2.oidx][info2.i];
    MQObject o0 = rv0.o;
    MQObject o1 = rv1.o;
    MQObject o2 = rv2.o;
    idx[0] = rv0.vidx;
    if(o0!=o2)
    {
      idx[2] = o0->AddVertex(rv2.p);
    } else idx[2] = rv2.vidx;
    if(o0!=o1)
    {
      idx[1] = o0->AddVertex(rv1.p);
    } else idx[1] = rv1.vidx;
    int vc = o0->GetVertexCount();
    int nfi = o0->AddFace(3, idx);
    if(nfi>=0)doc->AddSelectFace(info0.oidx, nfi);
  }
  return FALSE;
}
#endif //#ifdef ENABLE_TENUCHI

int FindVIndex(std::vector<int> &vidx, int vi, int start = 0)
{
  int num = vidx.size();
  for(int i=start;i<num;i++)
  {
    if(vidx[i]==vi)return i;
  }
  return -1;
}

int FindVIndexReverse(std::vector<int> &vidx, int vi, int start)
{
  for(int i=start;i>=0;i--)
  {
    if(vidx[i]==vi)return i;
  }
  return -1;
}
int FindVIndexReverse(std::vector<int> &vidx, int vi)
{
  return FindVIndexReverse(vidx, vi, vidx.size()-1);
}

bool SplitVertex(MQObject mqo, int target_fi, int del_vi, std::vector<int> &newtri_vidx, std::vector<MQCoordinate> &newtri_coord, std::vector<DWORD> &newtri_colors, int &matid)
{
  int numFV = mqo->GetFacePointCount(target_fi);
  if(numFV<=2)return false;
  
  std::vector<int> vidx(numFV);
  std::vector<MQPoint> pts(numFV);
  std::vector<MQCoordinate> coord(numFV);
  
  GetPointAndCoord(mqo, target_fi, numFV, vidx, pts, coord);

  std::vector<DWORD> colors(numFV);
  GetVertexColors(mqo, target_fi, numFV, colors);
  matid = mqo->GetFaceMaterial(target_fi);
  if(numFV==3)
  {
    int vi = FindVIndex(vidx, del_vi);
    int shift = 0;
    switch(vi)
    {
    case 0:
      shift = 2;
      break;
    case 1:
      shift = 0;
      break;
    case 2:
      shift = 1;
      break;
    default:
      return false;
    }
    std::rotate(vidx.begin(),vidx.begin()+shift,vidx.end());
    std::rotate(coord.begin(),coord.begin()+shift,coord.end());
    std::rotate(colors.begin(),colors.begin()+shift,colors.end());
    
    newtri_vidx = vidx;
    newtri_coord = coord;
    newtri_colors = colors;
    mqo->DeleteFace(target_fi);
    return true;
  } else {
    int vi = FindVIndex(vidx, del_vi);
    if(vi==-1)return false;
    
    int vi2 = FindVIndex(vidx, del_vi, vi+1);
    if(vi2!=-1)return false;
    
    int idx = vi-1;
    if(idx<0)idx = vidx.size()-1;
    std::rotate(vidx.begin(),vidx.begin()+idx,vidx.end());
    std::rotate(coord.begin(),coord.begin()+idx,coord.end());
    std::rotate(colors.begin(),colors.begin()+idx,colors.end());
    newtri_vidx = vidx;
    newtri_coord = coord;
    newtri_colors = colors;
    newtri_vidx.resize(3);
    newtri_coord.resize(3);
    newtri_colors.resize(3);
    
    vidx.erase(vidx.begin()+1);
    coord.erase(coord.begin()+1);
    colors.erase(colors.begin()+1);
    int newfi = mqo->AddFace(vidx.size(), &(vidx[0]));
    if(matid!=-1)mqo->SetFaceMaterial(newfi, matid);
    mqo->SetFaceCoordinateArray(newfi, &(coord[0]));
    SetVertexColors(mqo, newfi, numFV-1, colors);
    mqo->DeleteFace(target_fi);
    return true;
  }
}

void RemoveSpikeEdge(int del_vi, std::vector<int> &vidx, std::vector<MQCoordinate> &coord, std::vector<DWORD> &colors)
{
  int findstart = 0;
  int vi = -1;
  int numV = vidx.size();
  while(1)
  {
    if(vidx.size()<=3)break;
    vi = FindVIndex(vidx, del_vi, findstart);
    if(vi==-1)return;
    int li = LeftIndex(vi, numV);
    int ri = RightIndex(vi, numV);
    findstart = vi+1;
    if(vidx[li]==vidx[ri])
    {
      if(ri==vi+1)
      {
        vidx.erase(vidx.begin()+vi, vidx.begin()+(vi+2));
        coord.erase(coord.begin()+vi, coord.begin()+(vi+2));
        colors.erase(colors.begin()+vi, colors.begin()+(vi+2));
      } else {
        vidx.erase(vidx.begin()+(vi-1), vidx.begin()+vi+1);
        coord.erase(coord.begin()+(vi-1), coord.begin()+vi+1);
        colors.erase(colors.begin()+(vi-1), colors.begin()+vi+1);
      }
      findstart = vi;
    }
  }
}


template <typename VecType>
void EraseVector(std::vector<VecType>& vec, int start, int num)
{
  int vecSize = vec.size();
  if(num>vecSize || vecSize<=0 || num<=0)
  {
    vec.clear();
    return;
  }
  int end = start+num-1;
  if(end>=vecSize)end = vecSize-1;
  int eraseNum = end-start+1;
  if(eraseNum==1)vec.erase(vec.begin()+start);
  else           vec.erase(vec.begin()+start, vec.begin()+end+1);
  int remain = num - eraseNum;
  if(remain>0)
  {
    end = remain-1;
    vec.erase(vec.begin(), vec.begin()+end+1);
  }
}


/*

...
A      >>>>>>>
               del_vi
A      <<<<<<<
...

↓

...
A
...

*/
int __RemoveSpikeEdge2_TopDel(int del_vi, std::vector<int> &vidx, std::vector<MQCoordinate> &coord, std::vector<DWORD> &colors, int findstart, int numV)
{
  int li =  LeftIndex(findstart, numV);
  int ri = RightIndex(findstart, numV);
  const int l = vidx[li];
  const int r = vidx[ri];
  if(l==r)
  {
    EraseVector(vidx, li, 2);
    EraseVector(coord, li, 2);
    EraseVector(colors, li, 2);
    return findstart-1;//li;
  }
  return findstart+1;
}


/*

...
del_vi >>>>>>>
                B >>>>>>>
                          A
                B <<<<<<<
del_vi <<<<<<<
...

↓

...
...

*/
int __RemoveSpikeEdge2_NotTopDel(int del_vi, std::vector<int> &vidx, std::vector<MQCoordinate> &coord, std::vector<DWORD> &colors, int findstart, int numV, int maxComp)
{
  int li = findstart;
  int ri = findstart;
  for(int i=0;i<maxComp;i++)
  {
    li =  LeftIndex(li, numV);
    ri = RightIndex(ri, numV);
    const int l = vidx[li];
    const int r = vidx[ri];
    if(l!=r)break;
    if(l==del_vi)
    {
      i++;
      int numDel = i*2+1;
      if(i<maxComp && vidx[LeftIndex(li, numV)]==vidx[RightIndex(ri, numV)])numDel++;
      EraseVector(vidx, li, numDel);
      EraseVector(coord, li, numDel);
      EraseVector(colors, li, numDel);
      return findstart-i-1;//li-1;
    }
  }
  return findstart+1;
}

void RemoveSpikeEdge2(int del_vi, std::vector<int> &vidx, std::vector<MQCoordinate> &coord, std::vector<DWORD> &colors)
{
  int findstart = 0;
  int vi = -1;
  
  while(1)
  {
    int numV = vidx.size();
    if(numV<3 || findstart>=numV)return;
    int maxComp = (numV-1)/2;
    
    bool bTopDel = false;
    if(vidx[findstart]==del_vi)bTopDel = true;
    
    if(bTopDel) findstart = __RemoveSpikeEdge2_TopDel(del_vi, vidx, coord, colors, findstart, numV);
    else        findstart = __RemoveSpikeEdge2_NotTopDel(del_vi, vidx, coord, colors, findstart, numV, maxComp);
    if(findstart<0)findstart=0;
  }
}

int MakeFace2(MQObject mqoOut, std::vector<int> &vidx, std::vector<MQCoordinate> &coord, std::vector<DWORD> &colors, int matid)
{
  int addfi = mqoOut->AddFace(vidx.size(), &(vidx[0]));
  mqoOut->SetFaceCoordinateArray(addfi, &(coord[0]));
  SetVertexColors(mqoOut, addfi, colors.size(), colors);
  if(matid!=-1)mqoOut->SetFaceMaterial(addfi, matid);
  
  return addfi;
}

void _SplitVertex2_AddEdge()
{
}

bool SplitVertex2(MQObject mqo, int target_fi, int del_vi, MyGraph &graph, std::vector<int> &delFaceIdx)
{
  delFaceIdx.push_back(target_fi);//mqo->DeleteFace(target_fi);
  int numFV = mqo->GetFacePointCount(target_fi);
  if(numFV<=2)return false;
  
  std::vector<int> vidx(numFV);
  std::vector<MQPoint> pts(numFV);
  std::vector<MQCoordinate> coord(numFV);
  
  GetPointAndCoord(mqo, target_fi, numFV, vidx, pts, coord);

  std::vector<DWORD> colors(numFV);
  GetVertexColors(mqo, target_fi, numFV, colors);
  int matid = mqo->GetFaceMaterial(target_fi);
  if(numFV < 3)return false;
  
  RemoveSpikeEdge2(del_vi, vidx, coord, colors);
  
  int vi = FindVIndexReverse(vidx, del_vi);
  if(vi==-1)
  {
    MakeFace2(mqo, vidx, coord, colors, matid);
    return false;
  }
  
  //Ribbon Splitter >>>
  // del_vi == 3
  // vidx:  1 2 3 4 5 6 3 7 8 3 9
  //result >>>
  // edge1: 6 3 4 5
  // edge2: 8 3 7
  // edge3: 2 3 9 1
  while(1)
  {
    int vi_left = FindVIndexReverse(vidx, del_vi, vi-1);
    if(vi_left!=-1)
    {
      MyEdge *e = graph.new_edge();
      if(e!=NULL)
      {
        e->Set2(vidx, coord, colors, matid, vi_left, vi - vi_left);
        e->Rotate(e->v.size()-1);
      }
      vi = vi_left;
      EraseVector(vidx, vi_left, vi-vi_left);
    } else break;
  }
  
  int idx = LeftIndex(vi, vidx.size());
  std::rotate(vidx.begin(),vidx.begin()+idx,vidx.end());
  std::rotate(coord.begin(),coord.begin()+idx,coord.end());
  std::rotate(colors.begin(),colors.begin()+idx,colors.end());
  
  MyEdge *e = graph.new_edge();
  if(e!=NULL)e->Set(vidx, coord, colors, matid);
  return true;
}

int FindLeftVidx(std::vector<std::vector<int>> new_vidx, int find_vi)
{
  int n = new_vidx.size();
  
  for(int i=0;i<n;i++)
  {
    if(new_vidx[i][0]==find_vi)return i;
  }
  return -1;
}

int _CountMatid(std::vector<int> matids, int start)
{
  int num = matids.size();
  int ret = 1;
  int target = matids[start];
  for(int i=start+1;i<num;i++)
  {
    if(target!=matids[i])return ret;
    ret++;
  }
  return ret;
}

int GetMatId(std::vector<int> src_matids)
{
  int matid = -1;
  std::sort(src_matids.begin(), src_matids.end());
  int maxcnt = 0;
  int numMatId = src_matids.size();
  for(int i=0;i<numMatId;i++)
  {
    int cnt = _CountMatid(src_matids, i);
    int matid2 = src_matids[i];
    i+=cnt-1;
    if(cnt>maxcnt && matid2 != -1)
    {
      matid = matid2;
      maxcnt = cnt;
    }
  }
  
  return matid;
}

bool _MergeFace_Ring(MQObject mqoOut, std::vector<std::vector<int>> &src_vidx, std::vector<std::vector<MQCoordinate>> &src_coord, std::vector<std::vector<DWORD>> &src_colors, std::vector<int> &src_matids, int start = 0)
{
  int numRelFace = src_vidx.size();
  if(numRelFace<2)return false;
  
  std::vector<int> vidx;
  std::vector<MQCoordinate> coord;
  std::vector<DWORD> colors;
  std::vector<int> matids;
  
  vidx.push_back(src_vidx[start][0]);
  coord.push_back(src_coord[start][0]);
  colors.push_back(src_colors[start][0]);
  matids.push_back(src_matids[start]);
  int prev_end = src_vidx[start][2];
  src_vidx.erase(src_vidx.begin()+start);
  src_coord.erase(src_coord.begin()+start);
  src_colors.erase(src_colors.begin()+start);
  src_matids.erase(src_matids.begin()+start);
  
MergeFace_lp_start:
  if(src_vidx.size()>0)
  {
    int n = src_vidx.size();
    for(int i=0;i<n;i++)
    {
      if(prev_end==src_vidx[i][0])
      {
        prev_end = src_vidx[i][2];
        vidx.push_back(src_vidx[i][0]);
        coord.push_back(src_coord[i][0]);
        colors.push_back(src_colors[i][0]);
        matids.push_back(src_matids[i]);
        src_vidx.erase(src_vidx.begin()+i);
        src_coord.erase(src_coord.begin()+i);
        src_colors.erase(src_colors.begin()+i);
        src_matids.erase(src_matids.begin()+i);
        goto MergeFace_lp_start;
      }
    }
    if(prev_end!=vidx[0])return false;
  }
  if(vidx.size()<3 || prev_end!=vidx[0])return false;
  
  std::reverse(vidx.begin(),vidx.end());
  std::reverse(coord.begin(),coord.end());
  std::reverse(colors.begin(),colors.end());
  int addfi = mqoOut->AddFace(vidx.size(), &(vidx[0]));
  int new_matid = GetMatId(matids);
  if(new_matid!=-1)mqoOut->SetFaceMaterial(addfi, new_matid);
  mqoOut->SetFaceCoordinateArray(addfi, &(coord[0]));
  SetVertexColors(mqoOut, addfi, colors.size(), colors);
  return true;
}
bool _MergeFace_Arc(MQObject mqoOut, std::vector<std::vector<int>> &src_vidx, std::vector<std::vector<MQCoordinate>> &src_coord, std::vector<std::vector<DWORD>> &src_colors, std::vector<int> &src_matids, int start = 0)
{
  int numRelFace = src_vidx.size();
  if(numRelFace<2)return false;
  
  std::vector<int> vidx;
  std::vector<MQCoordinate> coord;
  std::vector<DWORD> colors;
  std::vector<int> matids;
  
  vidx.push_back(src_vidx[start][0]);
  coord.push_back(src_coord[start][0]);
  colors.push_back(src_colors[start][0]);
  matids.push_back(src_matids[start]);
  int prev_end = src_vidx[start][2];
  src_vidx.erase(src_vidx.begin()+start);
  src_coord.erase(src_coord.begin()+start);
  src_colors.erase(src_colors.begin()+start);
  src_matids.erase(src_matids.begin()+start);
  
MergeFace_lp_start:
  if(src_vidx.size()>0)
  {
    int n = src_vidx.size();
    for(int i=0;i<n;i++)
    {
      if(prev_end==src_vidx[i][0])
      {
        prev_end = src_vidx[i][2];
        vidx.push_back(src_vidx[i][0]);
        coord.push_back(src_coord[i][0]);
        colors.push_back(src_colors[i][0]);
        matids.push_back(src_matids[i]);
        src_vidx.erase(src_vidx.begin()+i);
        src_coord.erase(src_coord.begin()+i);
        src_colors.erase(src_colors.begin()+i);
        src_matids.erase(src_matids.begin()+i);
        goto MergeFace_lp_start;
      }
    }
    //if(prev_end!=vidx[0])return false;
  }
  if(vidx.size()<3/* || prev_end!=vidx[0]*/)return false;
  
  std::reverse(vidx.begin(),vidx.end());
  std::reverse(coord.begin(),coord.end());
  std::reverse(colors.begin(),colors.end());
  int addfi = mqoOut->AddFace(vidx.size(), &(vidx[0]));
  int new_matid = GetMatId(matids);
  if(new_matid!=-1)mqoOut->SetFaceMaterial(addfi, new_matid);
  mqoOut->SetFaceCoordinateArray(addfi, &(coord[0]));
  SetVertexColors(mqoOut, addfi, colors.size(), colors);
  return true;
}

int _FindRight(int vidxFind, std::vector<std::vector<int>> &src_vidx, int start = 0)
{
  int iFind = src_vidx[vidxFind][2];
  int num = src_vidx.size();
  for(int i=start;i<num;i++)
  {
    if(i==vidxFind)continue;
    if(src_vidx[i][0]==iFind)return i;
  }
  return -1;
}
int _FindLeft(int vidxFind, std::vector<std::vector<int>> &src_vidx, int start = 0)
{
  int iFind = src_vidx[vidxFind][0];
  int num = src_vidx.size();
  for(int i=start;i<num;i++)
  {
    if(i==vidxFind)continue;
    if(src_vidx[i][2]==iFind)return i;
  }
  return -1;
}

void _DeleteNotConnections(std::vector<std::vector<int>> &src_vidx, std::vector<std::vector<MQCoordinate>> &src_coord, std::vector<std::vector<DWORD>> &src_colors, std::vector<int> &src_matids, int start = 0)
{
  int num = src_vidx.size();
  for(int i=start;i<num;i++)
  {
    if(src_vidx[i][0]==src_vidx[i][2] || (_FindRight(i, src_vidx)==-1 && _FindLeft(i, src_vidx)==-1))
    {
      src_vidx.erase(src_vidx.begin()+i);
      src_coord.erase(src_coord.begin()+i);
      src_colors.erase(src_colors.begin()+i);
      src_matids.erase(src_matids.begin()+i);
    }
  }
}
void MergeFace(MQObject mqoOut, std::vector<std::vector<int>> src_vidx, std::vector<std::vector<MQCoordinate>> src_coord, std::vector<std::vector<DWORD>> src_colors, std::vector<int> src_matids)
{
  bool bRet = true;
  _DeleteNotConnections(src_vidx, src_coord, src_colors, src_matids);
  
  std::vector<std::vector<int>> vidx = src_vidx;
  std::vector<std::vector<MQCoordinate>> coord = src_coord;
  std::vector<std::vector<DWORD>> colors = src_colors;
  std::vector<int> matids = src_matids;
  
  do
  {
    vidx = src_vidx;
    coord = src_coord;
    colors = src_colors;
    matids = src_matids;
    bRet = _MergeFace_Ring(mqoOut, src_vidx, src_coord, src_colors, src_matids);
    if(bRet)
    {
      src_vidx = vidx;
      src_coord = coord;
      src_colors = colors;
      src_matids = matids;
    }
  }while(bRet);
  do
  {
    vidx = src_vidx;
    coord = src_coord;
    colors = src_colors;
    matids = src_matids;
    bRet = _MergeFace_Arc(mqoOut, src_vidx, src_coord, src_colors, src_matids);
    if(bRet)
    {
      src_vidx = vidx;
      src_coord = coord;
      src_colors = colors;
      src_matids = matids;
    }
  }while(bRet);
}
void UpdateGVIndex(std::vector<MyEdge> edges, std::vector<int> &table)
{
  int num = edges.size();
  for(int i=0;i<num;i++)
  {
    edges[i].UpdateGVIndex(table);
  }
}



void MergeFace2(MQObject mqoOut, MyGraph graph)
{
  graph.Prepara();
  
  std::vector<MyGraph> splitGraphByCon1;
  graph.SplitByConnection(splitGraphByCon1);
  
  std::vector<MyGraph> splitGraphByConMat1;
  int numGraph1 = splitGraphByCon1.size();
  for(int i=0;i<numGraph1;i++)
  {
    MyGraph &g = splitGraphByCon1[i];
    std::vector<MyGraph> tmpGraphs;
    g.SplitByMatid(tmpGraphs);
    splitGraphByConMat1.insert(splitGraphByConMat1.end(), tmpGraphs.begin(), tmpGraphs.end());
  }
  
  int N3 = splitGraphByConMat1.size();
  for(int i=0;i<N3;i++)
  {
    MyGraph &G = splitGraphByConMat1[i];
    G.MergeFace2_Ring(mqoOut);
  }
  
  //int N3 = splitGraphByConMat1.size();
  for(int i=0;i<N3;i++)
  {
    MyGraph &G = splitGraphByConMat1[i];
    G.MergeFace2_Arc(mqoOut);
  }
}

bool DeleteVertex(MQDocument doc, int del_oi, int del_vi)
{
  MQObject o = doc->GetObject(del_oi);
  if(o==NULL || o->GetLocking() || o->GetVisible()==0)return false;
  
  int numRelFace = o->GetVertexRelatedFaces(del_vi, NULL);
  if(numRelFace<=0)return false;
  
  std::vector<int> relFi;
  relFi.resize(numRelFace);
  o->GetVertexRelatedFaces(del_vi, &(relFi[0]));
  
  if(numRelFace==1)
  {
    o->DeleteFace(relFi[0]);
    return true;
  }
  
  std::vector<std::vector<int>> new_vidx(numRelFace);
  std::vector<std::vector<MQCoordinate>> new_coord(numRelFace);
  std::vector<std::vector<DWORD>> new_colors(numRelFace);
  std::vector<int> new_matids(numRelFace, -1);
  
  for(int i=0;i<numRelFace;i++)
  {
    std::vector<int> &newtri_vidx = new_vidx[i];
    std::vector<MQCoordinate> &newtri_coord = new_coord[i];
    std::vector<DWORD> &newtri_colors = new_colors[i];
    int &matid = new_matids[i];
    SplitVertex(o, relFi[i], del_vi, newtri_vidx, newtri_coord, newtri_colors, matid);
  }
  MergeFace(o, new_vidx, new_coord, new_colors, new_matids);
}

bool DeleteVertex2(MQDocument doc, int del_oi, int del_vi, bool bArcSplit)
{
  MQObject o = doc->GetObject(del_oi);
  if(o==NULL || o->GetLocking() || o->GetVisible()==0)return false;
  
  int numRelFace = o->GetVertexRelatedFaces(del_vi, NULL);
  if(numRelFace<=0)return false;
  
  std::vector<int> relFi;
  relFi.resize(numRelFace);
  o->GetVertexRelatedFaces(del_vi, &(relFi[0])); //同じ面番号が複数返る場合があるので注意！
  std::sort(relFi.begin(), relFi.end());
  relFi.erase(std::unique(relFi.begin(), relFi.end()), relFi.end());
  numRelFace = relFi.size();
  
  /*if(numRelFace==1)
  {
    o->DeleteFace(relFi[0]);
    return true;
  }*/
  
  //std::vector<MyEdge> new_edges(numRelFace);
  MyGraph new_edges;
  new_edges.SetArcSplit(bArcSplit);
  std::vector<int> delFaceIndex;
  
  for(int i=0;i<numRelFace;i++)
  {
    SplitVertex2(o, relFi[i], del_vi, new_edges, delFaceIndex);
  }
  if(new_edges.edges.size()>0)MergeFace2(o, new_edges);

  int numDelFace = delFaceIndex.size();
  for(int i=0;i<numDelFace;i++)
  {
    o->DeleteFace(delFaceIndex[i]);
  }
  return true;
}

BOOL TenTenPlugin::OnLeftButtonDown(MQDocument doc, MQScene scene, MOUSE_BUTTON_STATE& state)
{
  int mode = GetMode();
  HIT_TEST_PARAM param;
  switch(mode)
  {
#ifdef ENABLE_TENUCHI
  case 0:
  case 2:
    param.TestVertex = param.TestLine = FALSE;
    param.TestFace = TRUE;
    param.DisableFrontOnly = FALSE;
    HitTest(scene, state.MousePos, param);
    if(param.HitType==HIT_TYPE_FACE)
    {
      Triangulate1Poly(doc, param.ObjectIndex, param.FaceIndex);
      HitTest(scene, state.MousePos, param);
      
      BOOL bRet = AddVertexOnTri(doc, param.ObjectIndex, param.FaceIndex, param.HitPos, mode==0);
      if(mode==0 || !bRet)
      {
        if(mode!=0)
        {
      OutputDebugStringA("Undo3\n");
          RedrawAllScene();
          UpdateUndo();
        }
        return TRUE;
      }
      
      TriangulateSelected(doc);
      if(IsFrontTriSelected(doc, scene)==false)
      {
        errMsg = 1;
      OutputDebugStringA("Undo2\n");
        RedrawAllScene();
        UpdateUndo();
        return TRUE;
      }
      
      bRet = ReEdge(doc, scene);
      
      OutputDebugStringA("Undo\n");
      RedrawAllScene();
      UpdateUndo();
      OutputDebugStringA("UndoE\n");
      
      return TRUE;
    }
    break;
#endif //#ifdef ENABLE_TENUCHI
    
  case 1:
    param.TestFace = param.TestLine = FALSE;
    param.TestVertex = TRUE;
    param.DisableFrontOnly = FALSE;
    HitTest(scene, state.MousePos, param);
    if(param.HitType==HIT_TYPE_VERTEX)
    {
      bool bArcSplit = GetEnableArcSplit();
      DeleteVertex2(doc, param.ObjectIndex, param.VertexIndex, bArcSplit);
      
      RedrawAllScene();
      UpdateUndo();
      
      return TRUE;
    }
    break;
  default:
    return FALSE;
  }
  
  return FALSE;
}



