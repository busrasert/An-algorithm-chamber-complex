import sage.all
from sage.all import Polyhedron, matrix, identity_matrix
from sage.all import *

def createPolytope(A,b):
    Ab=[ [b[i]] + [ -A[i][j] for j in range(len(A[i])) ] for i in range(len(A))]
    P=Polyhedron(ieqs=Ab)
    return P

def intersection(x,y):
    return x.intersection(y)

def feasibility_cones(A,debug=0):
    """takes the matrix A and returns the set of cones obtained by fixing d hyperplanes """
    m=len(A)
    d=len(A[0])
    #if matrix(A).rank()==d:
    cones=[]
    for k in range (d,d+1):
        K=list(Subsets(range(m),m-k))#for fixing d hyperplanes, we put slack variables for m-d of them.
        for ell in K:
            Q=Polyhedron(lines=matrix(A).columns(), rays=[[1 if i==j else 0 for i in range(m)] for j in ell])
             #gives a cone for each vertex
            if Q not in cones:
                cones.append(Q)
    return cones
    #else:
     #   raise Exception("matrix A is not full rank")

def Intersection_Cones(cones):
    """takes cones (the set of feasibility cones) and returns the set of reduced intersections of subsets"""
    SI=[]#set of cones obtained by intersection
    for s in range((len(cones)),0,-1):
        Intersect=[reduce(intersection, S) for S in list(Subsets(cones,s))]#reduced intersection of subsets of cones
        for i in Intersect:
            SI.append(i)
    IC=[]#set of intersection cones containing only nonrepeated elements form IS
    for c in SI:
        if c not in IC:
            IC.append(c)
    return IC

def Chambers(IC,A,debug=0):
    """ takes intersection cones, gives the chambers of the complex"""
    Complex=[]
    d=len(A[0])
    m=len(A)
    #fx=[[1,1,1],[-1,-1,-1]]
    #bx=[(d*(d+1)/2),-(d*(d+1)/2)]
    #Q=createPolytope(fx,bx)
    At=(matrix(A)).transpose()
    Ker=At.right_kernel()
    #print Ker.gens()
    P=Polyhedron(lines=Ker.gens())
    for c in IC:
        if len(c.rays_list())>0:#this is for not to take the cone that has only two lines
            cR=[]#set of rays that has dimension bigger than 0 of a cone
            for r in c.rays_list():
                if debug==1:
                    print ("Ray of the cone: ", r)
                CP=createPolytope(A,r)
                if debug==2:
                    print ("Polytope created by the ray: ", CP)
                if CP.dimension()>0:
                    cR.append(CP)
            if debug==3:
                print ("Set of rays that give polytopes with dim > 0: ", cR)
            if len(cR)==len(c.rays_list()):
                Complex.append(c.intersection(P))
    Bn=[(a,b) for a in Complex for b in Complex]
    for x in Bn:
        if x[1]!=x[0] and x[1].dimension()==x[0].dimension():
            if x[0].intersection(x[1])==x[0]:
                try:
                    Complex.remove(x[1])
                except ValueError:
                    pass
    if debug==4:
        for c in Complex:
            #Cc=c.intersection(P)
            print ("Chamber of the Complex: ", c) #Coneref[c]
            print ("Rays: ", c.rays_list())
            #print "Lines: ", Cc.lines_list()
            for r in c.rays_list():
                print(createPolytope(A,r))#.intersection(Q))#.show())
    return Complex
