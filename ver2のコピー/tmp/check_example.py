from itertools import combinations
from math import comb
P=19
pts=[(0,0),(1,0),(0,1),(1,1),(2,15),(15,4),(11,15),(12,16)]
rows=[(3,[2,1,1,1,1,1,1,0]),(6,[3,2,2,2,2,2,2,2]),(1,[1,1,0,0,0,0,0,0]),(2,[1,0,1,1,1,1,0,0]),(3,[2,0,1,1,1,1,1,1]),(2,[1,1,1,1,1,0,0,0]),(6,[2,3,2,2,2,2,2,2]),(6,[2,2,3,2,2,2,2,2])]
def system(d,ps,ms):
    mons=[(i,j) for i in range(d+1) for j in range(d+1-i)]
    A=[]
    for (x,y),m in zip(ps,ms):
        for a in range(m):
            for b in range(m-a):
                A.append([(comb(i,a)*comb(j,b)*pow(x,i-a,P)*pow(y,j-b,P))%P if i>=a and j>=b else 0 for i,j in mons])
    return mons,A
def rr(A,n):
    A=[r[:] for r in A]; piv=[]; k=0
    for c in range(n):
        loc=next((t for t in range(k,len(A)) if A[t][c]%P),None)
        if loc is None:continue
        A[k],A[loc]=A[loc],A[k]
        v=pow(A[k][c],-1,P); A[k]=[(v*z)%P for z in A[k]]
        for t in range(len(A)):
            if t!=k:
                v=A[t][c]; A[t]=[(x-v*y)%P for x,y in zip(A[t],A[k])]
        piv.append(c);k+=1
    return A,piv
for k,d in [(3,1),(6,2)]:
    for inds in combinations(range(8),k):
        mons,A=system(d,[pts[i] for i in inds],[1]*k)
        assert len(rr(A,len(mons))[1])==len(mons)
for i in range(8):
    mons,A=system(3,pts,[2 if j==i else 1 for j in range(8)])
    assert len(rr(A,len(mons))[1])==len(mons)
print('Definition 3.4: all conditions verified over F19')
forms={}
def trim(a):
    while a and a[-1]%P==0:a.pop()
    return a
def gcd(a,b):
    a=trim(a[:]);b=trim(b[:])
    while b:
        while len(a)>=len(b):
            shift=len(a)-len(b);v=a[-1]*pow(b[-1],-1,P)%P
            for i,c in enumerate(b):a[i+shift]=(a[i+shift]-v*c)%P
            trim(a)
        a,b=b,a
    return a
for label,(d,ms) in enumerate(rows,2):
    mons,A=system(d,pts,ms); A,piv=rr(A,len(mons))
    free=[j for j in range(len(mons)) if j not in piv];assert len(free)==1
    v=[0]*len(mons);v[free[0]]=1
    for k,j in enumerate(piv):v[j]=-A[k][free[0]]%P
    form=[0]*(ms[0]+1)
    for c,(a,b) in zip(v,mons):
        if a+b==ms[0]:form[a]=c
    forms[label]=form
    print(label,'tangent cone coefficients of x^a y^(m-a):',form)
for i in [2,3]:
    f=forms[i];assert len(gcd(f,[(a*f[a])%P for a in range(1,len(f))]))==1
    assert f[-1]!=0
for i,j in [(4,5),(4,6),(3,7),(6,9),(6,8)]:
    assert len(gcd(forms[i],forms[j]))==1
    assert forms[i][-1] or forms[j][-1]
print('Both transverse pairs and all five empty triple intersections verified')
