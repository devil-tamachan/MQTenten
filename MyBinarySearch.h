
#ifndef _TAMABINSEARCH_
#define _TAMABINSEARCH_


int MyBinarySearch(int x, std::vector<int> &table)
{
  int l  = 0;
  int r = table.size()-1;
  
  int *_table = &(table[0]);
  
  while(l<r)
  {
    int c = (l+r)/2;
    if(_table[c] < x)
      l = c+1;
    else
      r = c;
  }
  if(_table[l]==x)return l;
  return -1;
}

#endif