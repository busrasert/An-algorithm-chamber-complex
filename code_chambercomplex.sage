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
    if debug==1:
        print ("num of vars:", d)
        print ("num of ineqs:", m)
    if matrix(A).rank()==d:
        if debug==2:
                print ("Full rank. Rank=", d)
        cones=[]
        for k in range (d,d+1):
            K=list(Subsets(range(m),m-k))#for fixing d hyperplanes, we put slack variables for m-d of them.
            for ell in K:
                Q=Polyhedron(lines=matrix(A).columns(), rays=[[1 if i==j else 0 for i in range(m)] for j in ell])#gives a cone for each vertex
                if debug==3:
                    print ("Signs for fixed rows: ", ell)
                    print ("Cone for fixed rows: ", Q)
                if Q not in cones:
                    cones.append(Q)
        if debug==4:
            print ("Set of all cones for fixed rows: ", cones)
        return cones
    else:
        raise Exception("matrix A is not full rank")

def Intersection_Cones(cones,nomults=1,debug=0):
    """takes cones (the set of feasibility cones) and returns the set of reduced intersections of subsets"""
    SI=[]#set of cones obtained by intersection
    for s in range((len(cones)),0,-1):
        Intersect=[reduce(intersection, S) for S in list(Subsets(cones,s))]#reduced intersection of subsets of cones
        if debug==1:
            print ("Subset of cones: ", s)
            print ("Intersection of the cones: ", intersect)
        for i in Intersect:
            SI.append(i)
    if nomults!=1:
        if debug==2:
            print ("Set of cones with possible multible repitations: ", SI)
        return SI
    IC=[]#set of intersection cones containes:contains onlly nonrepeated elements form IS
    for c in SI:
        if c not in IC:
            IC.append(c)
    if  debug==3:
        print ("Set of cones without repetitions: ", IC)
    return IC

def Chambers(IC,A,debug=0):
    """ takes intersection cones, gives the chambers of the complex"""
    Complex=[]
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
                #P=Polyhedron(lines=matrix(A).columns(), rays=c.rays_list())
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
            Cc=c.intersection(P)
            print ("Chamber of the Complex: ", Cc) #,Coneref[c]
            print ("Rays: ", Cc.rays_list())
            #print ("Lines: ", Cc.lines_list())
            #for r in Cc.rays_list():
                #createPolytope(A,r).show()

    return Complex
A=[[1,0],[0,1],[-1,0],[0,-1],[1,1]]
#A=[[-1,0],[0,-1],[1,1]]
#A=[[-1,1],[1,1],[1,0],[2,-1],[-1,0]]
#A=[[0,-1],[-1,0],[4,11],[7,13],[3,5]]
#A=[[0,1,1],[1,0,1],[1,0,0],[0,1,0],[0,0,-1],[0,0,1],[0,-1,0],[1,-1,1],[-1,0,0],[-1,0,1]]
#A=[[0, 1, 0, 0, 1, 1, 0, 0, 0, 0], [1, 0, 1, 0, 0, 0, 1, 0, 0, 0], [0, 1, 0, 1, 0, 0, 0, 1, 0, 0], [0, 0, 1, 0, 1, 0, 0, 0, 1, 0], [1, 0, 0, 1, 0, 0, 0, 0, 0, 1], [1, 0, 0, 0, 0, 0, 0, 1, 1, 0], [0, 1, 0, 0, 0, 0, 0, 0, 1, 1], [0, 0, 1, 0, 0, 1, 0, 0, 0, 1], [0, 0, 0, 1, 0, 1, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0, 1, 1, 0, 0]]
cones=feasibility_cones(A)
IC=Intersection_Cones(cones,1,0)
Chambers(IC,A,4)
