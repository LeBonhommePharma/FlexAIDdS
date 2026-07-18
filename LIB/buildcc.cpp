#include "flexaid.h"
#include <math.h>
#include <cmath>
#include <limits>
/******************************************************************************
 * SUBROUTINE buildcc builds the cartesian coordinates of the tot atoms present
 * in array list according to the reconstruction data.
 *
 * Fail-closed: singular GPA frames or non-finite outputs return false and
 * write NaNs so callers cannot treat stale coordinates as a successful rebuild.
 ******************************************************************************/
bool buildcc(FA_Global* FA,atom* atoms,int tot,int list[]){
  int an,i,j;
  float x[4],y[4],z[4];
  float a,b,c,op,cx,cy,cz,d,xn,yn,zn,ct,st,xk,yk,zk,angPI,dihPI;
  bool ok = true;
  const float nanf = std::numeric_limits<float>::quiet_NaN();

  for(an=0;an<tot;an++){
    
    for(i=1;i<=3;i++)
    {
      j=atoms[list[an]].rec[i-1];
      if(j != 0)
      {
      	x[i]=atoms[j].coor[0];
      	y[i]=atoms[j].coor[1];
      	z[i]=atoms[j].coor[2];
      }
      else if(i==1)
      {
      	x[i]=1.0f+FA->ori[0];
      	y[i]=0.0f+FA->ori[1];
      	z[i]=0.0f+FA->ori[2];
      }
      else if(i==3)
      {
      	x[i]=0.0f+FA->ori[0];
      	y[i]=1.0f+FA->ori[1];
      	z[i]=0.0f+FA->ori[2];
      }
      else
      {
      	x[i]=0.0f+FA->ori[0];
      	y[i]=0.0f+FA->ori[1];
      	z[i]=0.0f+FA->ori[2];
      }	

      // perturb atom coordinates
      x[i]+=1e-10f;
      y[i]+=1e-10f;
      z[i]+=1e-10f;
    }
    
    a=y[1]*(z[2]-z[3])+y[2]*(z[3]-z[1])+y[3]*(z[1]-z[2]);
    b=z[1]*(x[2]-x[3])+z[2]*(x[3]-x[1])+z[3]*(x[1]-x[2]);
    c=x[1]*(y[2]-y[3])+x[2]*(y[3]-y[1])+x[3]*(y[1]-y[2]);
    op=sqrtf(a*a+b*b+c*c);
    // Collinear / singular GPA frame: fail closed (NaN + status).
    if (!(op > 1e-12f) || !std::isfinite(op)) {
      atoms[list[an]].coor[0] = nanf;
      atoms[list[an]].coor[1] = nanf;
      atoms[list[an]].coor[2] = nanf;
      ok = false;
      continue;
    }

    cx=a/op;
    cy=b/op;
    cz=c/op;

    a=x[2]-x[1];
    b=y[2]-y[1];
    c=z[2]-z[1];

    const float ref_len = sqrtf(a*a+b*b+c*c);
    if (!(ref_len > 1e-12f) || !std::isfinite(ref_len)) {
      atoms[list[an]].coor[0] = nanf;
      atoms[list[an]].coor[1] = nanf;
      atoms[list[an]].coor[2] = nanf;
      ok = false;
      continue;
    }
    d=1.0f/ref_len;
    op=atoms[list[an]].dis*d;
    xn=a*op;
    yn=b*op;
    zn=c*op;

    a=cx*cx;
    b=cy*cy;
    c=cz*cz;

    angPI = (float)(atoms[list[an]].ang*PI/180.0f);
    st = -sinf(angPI);
    ct = cosf(angPI);

    op=1.0f-ct;

    xk=(cx*cz*op-cy*st)*zn+((1.0f-a)*ct+a)*xn+(cx*cy*op+cz*st)*yn;
    yk=(cy*cx*op-cz*st)*xn+((1.0f-b)*ct+b)*yn+(cy*cz*op+cx*st)*zn;
    zk=(cz*cy*op-cx*st)*yn+((1.0f-c)*ct+c)*zn+(cz*cx*op+cy*st)*xn;
    
    dihPI = (float)(atoms[list[an]].dih*PI/180.0f);
    st = sinf(dihPI);
    ct = cosf(dihPI);

    op=1.0f-ct;

    cx=(x[2]-x[1])*d;
    cy=(y[2]-y[1])*d;
    cz=(z[2]-z[1])*d;   
    
    a=cx*cx;
    b=cy*cy;
    c=cz*cz;

    x[0]=(cx*cz*op-cy*st)*zk+((1.0f-a)*ct+a)*xk+(cx*cy*op+cz*st)*yk+x[1];
    y[0]=(cy*cx*op-cz*st)*xk+((1.0f-b)*ct+b)*yk+(cy*cz*op+cx*st)*zk+y[1];
    z[0]=(cz*cy*op-cx*st)*yk+((1.0f-c)*ct+c)*zk+(cz*cx*op+cy*st)*xk+z[1];

    if (!std::isfinite(x[0]) || !std::isfinite(y[0]) || !std::isfinite(z[0])) {
      atoms[list[an]].coor[0] = nanf;
      atoms[list[an]].coor[1] = nanf;
      atoms[list[an]].coor[2] = nanf;
      ok = false;
      continue;
    }
    
    atoms[list[an]].coor[0]=x[0];
    atoms[list[an]].coor[1]=y[0];
    atoms[list[an]].coor[2]=z[0];
   }

  return ok;
}
