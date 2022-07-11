
import numpy as np
import time
start_time = time.time()



f = open('./gcode/embedded hexagon_0.2mm_PLA_MK3S_7m.gcode','r') #reads in with type _io.TextIOWrapper which is a not super usable object

f = f.read() #The read function returns f as a str which has much more usability

#print(type(f)) #debug

#print(f) #debug
li = list(f.split("\n")) #\n is the carriage return of each line
#print(type(li))

perimeters = [] #initializes filtered as a list
 
for x in li: #iterates through the list storing each iteration of li into x 
    
    #print(x)
    if (x.find("; perimeter")!= -1 and x.find("; perimeters")== -1) or x.find("; restore layer Z")!= -1 or x.find("; move to first perimeter point") !=-1: #checks if the string x has a substring ; perimeter and also checks that ; perimeters is not included as well
        perimeters.append(x) #appends x to the filtered list
print(len(f))

print("filtered length:")
print(len(perimeters))



#
#iterates through perimeters to output
#

z = 0.0 #initializes z = 0 
ztest = []
xyz = [] #initializes the xyz 2D array 
new = [] #initializes the individual line

for string in perimeters:
    intermed1 = list(string.split(" "))
    #print(intermed1)
    for string2 in intermed1:
        #print(string2)
        if string2.find("Z") !=-1 and len(string2) != 1:
            #print(string2)
            
            z = float(string2.replace("Z",""))
            ztest.append(z)
        elif string2.find("X") !=-1:
            string3 = string2.replace("X","")
            #print(string3)
            new.append(float(string3))
            #print(new)
        elif string2.find("Y") !=-1:
            string3 = string2.replace("Y","")
            #print(string3)
            new.append(float(string3))
            new.append(z)
            #print(new)
            xyz.append(new)    

                
    new = []
#print(ztest)

xyz = np.array(xyz, dtype = 'float32')
#print(xyz[:,2])



#Variable declaration for hexagon finding
lastPt = np.array([],dtype = 'float32')
delta = np.array([],dtype = 'float32')
curmag = 0.0
lastmag = 0.0
lastdelta = np.array([],dtype = 'float32')
hexagons = np.array([],dtype = 'float32')
angle = 0.0
hexcount = 0
fullHexCount = 0
hmmArray = np.array([[0,0],[0,0]],dtype = 'float32')
hmmNewline = np.array([],dtype = 'float32')
hexagonCtr = 0.0
HexagonCtrs = np.array([0,0,0],dtype = 'float32')
HexagonDia = np.array([0],dtype = 'float32')
HexagonAngle = np.array([0],dtype = 'float32')
ctrfindingMatrix = np.empty((3,3),dtype = 'float32')
hexagonPoints = np.empty((3,5),dtype = 'float32')

def mag(delta):
    return np.sqrt(delta[0]**2+delta[1]**2+delta[2]**2)

def minor(Marray,column):
    front = Marray[:,:column-1]
    #print(front)
    back = Marray[:,column:]
    #print(back)
    det = np.append(front,back,axis=1)
    #print(det)
    return np.linalg.det(det)


def centerFind(hexagon):
    #print(hexagon)
    midpoints = hexagon[1:4,:]
    #print(midpoints)
    newline = np.empty((1,4),dtype='float32')
    centerMatrix = np.empty((3,4),dtype='float32')
    for i in range(3):
        
        #print("i")
        newline = [midpoints[i,0]**2+midpoints[i,1]**2,midpoints[i,0],midpoints[i,1],1]
        #print(newline)
        
        centerMatrix[i,:] = newline
        #print(centerMatrix)

    xc = .5 * minor(centerMatrix,2)/minor(centerMatrix,1)
    yc =-.5 * minor(centerMatrix,3)/minor(centerMatrix,1)
    return np.array([xc,yc,0])

for point in xyz:
    if lastPt.size >0:
        delta = point-lastPt
        curmag = mag(delta)
        hmmNewline=np.append(hmmNewline,round(curmag,1))
        if lastmag != 0.0:
            angle  = 180*np.arccos(np.dot(delta[:2],lastdelta[:2])/(curmag*lastmag))/np.pi
            hmmNewline=np.append(hmmNewline,angle)
            angle = round(angle,1)
            #print(angle)
            #print(abs(curmag-lastmag))
            if abs(angle-60) <.1:
                hexcount +=1
                hexagonPoints[:,hexcount-1] = point

            else:
                hexcount = 0
                hexagonPoints = np.empty((3,5))
            

            if hexcount == 4:
                fullHexCount +=1
                #print(lastdelta)
                #print(lastmag)
                #print([lastmag*np.cos(2*np.pi/3),lastmag*np.sin(2*np.pi/3),0])
                #.45*np.ones(np.size(HexagonDia),dtype = 'float32').45*np.ones(np.size(HexagonDia),dtype = 'float32')print(lastmag)
                HexagonAngle = np.vstack([HexagonAngle,np.mod(round(180/np.pi*np.arctan2(lastdelta[1],lastdelta[0]),1),60)])
                
                hexagonPoints = np.transpose(hexagonPoints)
                #print(hexagonPoints)
                hexagonCtr = centerFind(hexagonPoints)
                hexagonCtr[2] = point[2]
                
                HexagonCtrs = np.vstack([HexagonCtrs,hexagonCtr])
                hexcount = 0
                hexagonPoints = np.empty((3,5))
                HexagonDia = np.vstack([HexagonDia,round(lastmag*np.sqrt(3),3)])
            #print(angle)
            #print(fullHexCount)
        else:
            hmmNewline=np.append(hmmNewline,0)

        #print(hmmNewline)
        hmmArray=np.vstack([hmmArray,hmmNewline]) 
        hmmNewline = []

        #print(delta)
        lastdelta = delta
        lastmag = curmag
        
    
    lastPt = point
#print(hmmArray)        
#print(np.shape(hmmArray))
#print(fullHexCount)
#print(HexagonCtrs)

HexagonEquivalent = np.empty((np.size(HexagonDia),1),dtype = 'float32')
HexagonEquivalent = HexagonDia -.45*np.ones((np.size(HexagonDia),1),dtype = 'float32')
#print(HexagonEquivalent)
#print(HexagonAngle)
#print(np.shape(HexagonEquivalent))
#print(np.shape(HexagonCtrs))
#print(np.shape(HexagonEquivalent))


HexagonFull = np.append(HexagonCtrs,HexagonEquivalent, axis=1)
HexagonFull = np.append(HexagonFull,HexagonAngle, axis=1)
HexagonFull = HexagonFull[1:,:].copy()
#print(HexagonFull,flush=True)
#print(np.shape(HexagonFull))

#print("--- %s seconds ---" % (time.time() - start_time))


fullNuts = np.empty((2,5),dtype='float32')

for hexagonfrag in HexagonFull:
    noMatch = False
    index =0
    if np.linalg.norm(fullNuts) !=0:
        for i in range(len(fullNuts[:,0])):
                hexagonfragNoZ = np.append(hexagonfrag[:2],hexagonfrag[3:])
                CurrentEntry = fullNuts[i,:] 
                CurrentEntryNoZ = np.append(CurrentEntry[:2],CurrentEntry[3:])
                print(hexagonfragNoZ)
                #print(np.shape([hexagonfrag]))
                #print(np.linalg.norm(hexagonfragNoZ-lastFragNoZ))
                #print(np.linalg.norm(fullNuts[i,:]-hexagonfrag))
                if np.linalg.norm(hexagonfragNoZ-CurrentEntryNoZ) !=0:
                    noMatch = True
                    index = i
                    #print(np.shape(fullNuts[i,:]))
        if noMatch: 
            fullNuts[index,:] = hexagonfrag
        else:
            fullNuts = np.vstack((fullNuts,hexagonfrag))
        noMatch = False
    else:
        fullNuts = hexagonfrag
        fullNuts = fullNuts[0,:]

        #print(np.shape(fullNuts))
             

print(np.shape(fullNuts))
        
