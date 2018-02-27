def udiff(F,u,v):
  """
  This method calculates the partial derivative of function F with respect to u.
  """
  return (F(u+10**(-6),v)-F(u,v))/(10**(-6))

def vdiff(F,u,v):
  """
  This method calculates the partial derivative of function F with respect to v.
  """
  return (F(u,v+10**(-6))-F(u,v))/(10**(-6))

def f(G,X,Y,Z,a,b):
  """
  This method calculates the value of ||ruXrv|| at point (a,b).
  """
  xu=udiff(X,a,b)
  xv=vdiff(X,a,b)
  yu=udiff(Y,a,b)
  yv=vdiff(Y,a,b)
  zu=udiff(Z,a,b)
  zv=vdiff(Z,a,b)
  cx=((yu*zv-yv*zu)**2)
  cy=((xu*zv-xv*zu)**2)
  cz=((xu*yv-xv*yu)**2)
  return G(a,b)*((cx+cy+cz)**0.5)

def surfintegrate(G,X,Y,Z,U,V):
  """
  This method calulates the value of surface integral by the application of Simpson's rule for double integrals.
  """
  N=100
  M=100
  h=(U[1]-U[0])/(2*M)
  k=(V[1]-V[0])/(2*N)
  u=U[0]
  v=V[0]
  I=0
  I+=(f(G,X,Y,Z,U[0],V[0])+f(G,X,Y,Z,U[0],V[1])+f(G,X,Y,Z,U[1],V[0])+f(G,X,Y,Z,U[1],V[1]))*h*k/9
  j=1
  S=[]
  s=0
  while j<=N:
    s+=f(G,X,Y,Z,U[0],v+(2*j-1)*h)
    j+=1
  S.append(4*s)
  s=0
  j=1
  while j<=N-1:
    s+=f(G,X,Y,Z,U[0],v+(2*j)*h)
    j+=1
  S.append(2*s)
  s=0
  j=1
  while j<=N:
    s+=f(G,X,Y,Z,U[1],v+(2*j-1)*h)
    j+=1
  S.append(4*s)
  s=0
  j=1
  while j<=N-1:
    s+=f(G,X,Y,Z,U[1],v+(2*j)*h)
    j+=1
  S.append(2*s)
  s=0
  i=1
  while i<=M:
    s+=f(G,X,Y,Z,u+(2*i-1)*k,V[0])
    i+=1
  S.append(4*s)
  s=0
  i=1
  while i<=M-1:
    s+=f(G,X,Y,Z,u+2*i*k,V[0])
    i+=1
  S.append(2*s)
  s=0
  i=1
  while i<=M:
    s+=f(G,X,Y,Z,u+(2*i-1)*k,V[1])
    i+=1
  S.append(4*s)
  s=0
  i=1
  while i<=M-1:
    s+=f(G,X,Y,Z,u+2*i*k,V[1])
    i+=1
  S.append(2*s)
  s=0
  j=1
  while j<=N:
    i=1
    while i<=M:
      s+=f(G,X,Y,Z,u+(2*i-1)*k,v+(2*j-1)*h)
      i+=1
    j+=1
  S.append(16*s)
  s=0
  j=1
  while j<=N-1:
    i=1
    while i<=M:
      s+=f(G,X,Y,Z,u+(2*i-1)*k,v+2*j*h)
      i+=1
    j+=1
  S.append(8*s)
  s=0
  j=1
  while j<=N:
    i=1
    while i<=M-1:
      s+=f(G,X,Y,Z,u+2*i*k,v+(2*j-1)*h)
      i+=1
    j+=1
  S.append(8*s)
  s=0
  j=1
  while j<=N-1:
    i=1
    while i<=M-1:
      s+=f(G,X,Y,Z,u+2*i*k,v+2*j*h)
      i+=1
    j+=1
  S.append(4*s)
  return (sum(S)*h*k/9)+I

"""
Main Code
"""
import numpy
G=lambda u,v:1
X=lambda u,v:(2+numpy.cos(u))*numpy.cos(v)
Y=lambda u,v:(2+numpy.cos(u))*numpy.sin(v)
Z=lambda u,v:numpy.sin(u)
U=[0,2*numpy.pi]
V=[0,2*numpy.pi]
J=surfintegrate(G,X,Y,Z,U,V)
print(J)
