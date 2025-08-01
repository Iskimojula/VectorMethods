'''
name: VectorMethods.py
Author: Liu Yang
创建日期: 2025-07-10

Pseidocode
Algorithm: Vector Methods for calculating segmentation points {xij, kij}
Input: N, S, a, D

if i ==1: x = mod(vector: OS), k = 1
for i ∈ [2, N]
  coordinate: Pip = [sum(D[:i])*cos(π/6), sum(D[:i])*sin(π/6)]
  coordinate: Pin = [sum(D[:i])*cos(π/6), -sum(D[:i])*sin(π/6)]
  vector: PipPin = sub(vector: OPin, vector: OPip)
  deta = mod(vector: PipPin)/(i-1)

  for j in range(i -1):
vector: PipPij = dotproduct(vector: PipPin, deta*j/mode(vector: PipPin))
coordinate: Pij = sum(vector: OPip, vector: PipPij)
if j == 0: vector: QPij0 = [3*a/4, sqrt(3)*a/4]; else : QPij0 = [a,0]
for quandrant ∈ [0,5]:
  rotateangle = quandrant*π/3
  vector: OPijq = rotate(vector: OPij, rotateangel)
  vector: QPijq = rotate(vector: QPij0, rotateangel)
  vector: SQijq = vector: OPijq - vector: QPijq – vector: OS
  xijq = mode(vector: SQijq)
  kijq = 1

xijq ∈ X and sort(X)
xij, Δkij = unique(X) and return count
kij = sum(Δkij[:i+1])c

return xij, kij

'''


import numpy
import math
import random

#input unit:um

# a is the length of hexagon pixel side
# Di is the distance between adjacent pixels
# I is the number of pixels on axis
I = 12
a = 0 #600/math.sqrt(3)
D = [0,650,650,650,650,650,650,650,650,650,650]

theta = math.radians(30)

# S is center of contact process
def genS(D_ave):
    A=[0,0]
    B=[D_ave*math.cos(theta),D_ave*math.sin(theta)]
    C=[D_ave*math.cos(theta),0]
    u = random.random()
    v = random.random()
    if u + v > 1:
        u, v = 1 - u, 1 - v
    w = 1 - u - v
    S = (
        u * A[0] + v * B[0] + w * C[0],
        u * A[1] + v * B[1] + w * C[1]
    )
    return S
S = genS(sum(D)/(len(D)-1))
S = [0,0]

#Calculate the distance between two points.
def dist(v1,v2):
    return math.sqrt((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1]))

#Calculate the mod of the vector
def vecmod(v1):
    return math.sqrt(v1[0]*v1[0]+v1[1]*v1[1])

#Calculate the sum of v1 and v2
def vecsum(v1,v2):
    return [v1[0]+v2[0],v1[1]+v2[1]]

#Calculate the sub of v1 and v2
def vecsub(v1,v2):
    return [v1[0]-v2[0],v1[1]-v2[1]]

#Calculate the dotproduct of v1 and valuec
def vecdotprocduct(v1,value):
    return [v1[0]*value,v1[1]*value]

#Rotate the v1 counterclockwise by theta degrees
def vecrotate(v1,theta): ##deg
    return [(math.cos(theta)*v1[0] - math.sin(theta)*v1[1]),(math.sin(theta)*v1[0] + math.cos(theta)*v1[1])]

def printlist(list):
    for l in list:
        print(l)

#Calculate vector:OPij, QPij in different quadrant
class Hexagon(object):
    def __init__(self,P,label):
        self.P0 = P
        self.label = label
        self.Plist = []
        self.Plist.append(P)
        self.Plist.append(vecrotate(P,math.radians(60)))
        self.Plist.append(vecrotate(P,math.radians(120)))
        self.Plist.append(vecrotate(P,math.radians(180)))
        self.Plist.append(vecrotate(P,math.radians(240)))
        self.Plist.append(vecrotate(P,math.radians(300)))
        self.K = len(self.Plist)

        self.QPlist = []
        
       
        if label == 'in_axis':
            QP0 = [3*a/4,math.sqrt(3)*a/4]
        else:
            QP0 = [a,0]
        self.QPlist.append(QP0)
        self.QPlist.append(vecrotate(QP0,math.radians(60)))
        self.QPlist.append(vecrotate(QP0,math.radians(120)))
        self.QPlist.append(vecrotate(QP0,math.radians(180)))
        self.QPlist.append(vecrotate(QP0,math.radians(240)))
        self.QPlist.append(vecrotate(QP0,math.radians(300)))
        

        self.R = vecmod(vecsub(self.P0,self.QPlist[0]))
        self.Qlist = []
        self.SQlist = []
        self.Qmapxyzlist=[] #Qx,  Qy, mod(SQ)
        for cnt in range(self.K):
            self.Qlist.append(vecsub(self.Plist[cnt],self.QPlist[cnt]))
            self.SQlist.append(vecsub(vecsub(self.Plist[cnt],S),self.QPlist[cnt]))
            Qmap = [0,0,0]
            Qmap[0] = vecsub(self.Plist[cnt],self.QPlist[cnt])[0]
            Qmap[1] = vecsub(self.Plist[cnt],self.QPlist[cnt])[1]
            Qmap[2] = vecmod(vecsub(vecsub(self.Plist[cnt],S),self.QPlist[cnt]))
            self.Qmapxyzlist.append(Qmap)
    
    def print_Q(self):
        for Q in self.Qlist:
            print(Q[0],Q[1])
    def print_P(self):
        for P in self.Plist:
            print(P[0],P[1])
    def print_Qmap(self):
        for Qmap in self.Qmapxyzlist:
            print(Qmap[0],Qmap[1],Qmap[2])


#i = 1
x_list = []
P0 = [0,0]
Pinit = Hexagon(P0,'in_axis')
#pixel_dict = {'1':Pinit} # for debug
x_list.append(vecmod(S))

i=2
while i <= I:
    #calculate coordinate: Pip,Pin and vector:PipPin
    rou = sum(D[:i])
    Pip = [rou*math.cos(theta),rou*math.sin(theta)]
    Pin = [rou*math.cos(theta),-rou*math.sin(theta)]
    vec_PipPin = vecsub(Pin,Pip)
    #list = [] # for debug
    deta = dist(Pip,Pin)/(i-1)

    #pixels in quadrant 0
    for k in range(i-1):
        vec_PipPik = vecdotprocduct(vec_PipPin,deta*k/vecmod(vec_PipPin))
        Pik = vecsum(Pip,vec_PipPik)
        if k == 0 : label = 'in_axis' 
        else: label = 'out_axis'

        #calculate vector:OPij, QPij in other quadrant
        pixel = Hexagon(Pik,label)
        #list.append(pixel) #for debug
        #pixel.print_P() #for debug
        #pixel.print_Qmap() #for debug
        for SQ in pixel.SQlist:
            x_list.append(vecmod(SQ))
            #print(vecmod(SQ)) #fordebug
            
    i+=1
    
#printlist(x_list)
x_list = sorted([round(x,5) for x in x_list])

r,kij = numpy.unique(numpy.array(x_list),return_counts=True)
k_list = []
for i in range(len(kij)):
    k_list.append(sum(kij[:i+1]))

#print result
if 1:
    print("x-k",len(r))
    for i in range(len(r)):
        print(r[i],k_list[i])
