
import numpy as np
import time
start_time = time.time()



f = open('./gcode/embedded hexagon_0.2mm_PLA_MK3S_6m.gcode','r') #reads in with type _io.TextIOWrapper which is a not super usable object

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

xyz = [] #initializes the xyz 2D array 
new = [] #initializes the individual line

for string in perimeters: #iterates through the perimeter points with string containing each point
    intermed1 = list(string.split(" ")) #string.split breaks a string into a list of substrings
    #print(intermed1)
    for string2 in intermed1:   #iterates through all substrings and reads each string looking for x,y, or z 
        #print(string2)
        if string2.find("Z") !=-1 and len(string2) != 1: #checks if the substring is like Z0.6 compared to Z
            #print(string2)
            
            z = float(string2.replace("Z","")) #removes Z from Z.6 to make .6
            
        elif string2.find("X") !=-1:    #checks for X 
            string3 = string2.replace("X","") #removes X to go from X254 to 254
            #print(string3)
            new.append(float(string3))  #adds the x value 
            #print(new)
        elif string2.find("Y") !=-1:    #checks for Y
            string3 = string2.replace("Y","") #removes Y to go from Y254 to 254
            #print(string3)
            new.append(float(string3)) #adds y to the new array
            new.append(z)               #adds z to create an array of x,y,z
            #print(new)
            xyz.append(new)     #adds x,y,z point to the xyz array

                
    new = [] #emptys new array
#print(ztest)

xyz = np.array(xyz, dtype = 'float32') #converts xyz into a numpy array
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

def mag(delta): #mag function call to evaluate the magnitude 
    return np.sqrt(delta[0]**2+delta[1]**2+delta[2]**2) #basic magnitude function

def minor(Marray,column): #minor function. Required for evaluating a circle center. Technically this is an incomplete minor function bc a 
                            #minor function is the determinant of an array with one row and one collumn removed
                            #in this case, the first row is always removed and the column arg is the collumn that is removed
    front = Marray[:,:column-1] #pulls front half of array
    #print(front)
    back = Marray[:,column:] #pulls back half of array
    #print(back)
    det = np.append(front,back,axis=1) #combines both vectors. could have used np.hstack instead of np.append
    #print(det)
    return np.linalg.det(det) #returns the determinant of the minored array. single value with type float


def centerFind(hexagon): #center finding function which calls minor function. Algorithm is outlined in https://math.stackexchange.com/questions/213658/get-the-equation-of-a-circle-when-given-3-points
    #print(hexagon)
    midpoints = hexagon[1:4,:]  #the first and last hexagon vector are both shorter than the rest of the vectors because of slicer stuff
    #print(midpoints)
    centerMatrix = np.empty((3,4),dtype='float32') #creates matrix outlined in the above link. 
    for i in range(3):
        
        #print("i")
        centerMatrix[i,:] = [midpoints[i,0]**2+midpoints[i,1]**2,midpoints[i,0],midpoints[i,1],1] # stores valuse as x^2 + y^2, x, y , 1
        #print(centerMatrix)

    xc = .5 * minor(centerMatrix,2)/minor(centerMatrix,1) #follows algorithm of the above link
    yc =-.5 * minor(centerMatrix,3)/minor(centerMatrix,1) # ^
    return np.array([xc,yc,0]) #returns a concatenated array of the x and y center points

for point in xyz: #iterates through xyz points and breaks it up with point
    if lastPt.size >0: #checks that it is the second iteration
        delta = point-lastPt #finds the delta vector between the last point and the current point
        curmag = mag(delta) #evaluates current vector magnitude
        hmmNewline=np.append(hmmNewline,round(curmag,1)) #adds the rounded magnitude of the current vector
        if lastmag != 0.0: #verifies that this is the 3rd iteration
            angle  = 180*np.arccos(np.dot(delta[:2],lastdelta[:2])/(curmag*lastmag))/np.pi #evaluates angle with SOH CAH TOA
            hmmNewline=np.append(hmmNewline,angle) #adds angle to the hmmnewline array
            angle = round(angle,1) #rounds off angle to check for a hexagon match
            #print(angle)
            #print(abs(curmag-lastmag))
            if abs(angle-60) <.1: #checks if the absolute error between the angle and 60 degrees is less than .1
                hexcount +=1 #hexcount increments
                hexagonPoints[:,hexcount-1] = point #stores the current point into a hexagonPoints array

            else:
                hexcount = 0 #restarts the hex count bc the error is too high
                hexagonPoints = np.empty((3,5)) #emptys the hexagon points array
            

            if hexcount == 4: #if the hexagon count is 4... it seems that there us an error with the angle on the first and the last vector so hexcount has a max val of 4 for a proper hexagon
                fullHexCount +=1 #adds 1 to full hex cound for debug purposes
                #print(lastdelta)
                #print(lastmag)
                #print([lastmag*np.cos(2*np.pi/3),lastmag*np.sin(2*np.pi/3),0])
                #.45*np.ones(np.size(HexagonDia),dtype = 'float32').45*np.ones(np.size(HexagonDia),dtype = 'float32')print(lastmag)
                HexagonAngle = np.vstack([HexagonAngle,np.mod(round(180/np.pi*np.arctan2(lastdelta[1],lastdelta[0]),1),60)]) #rounds angle to be less than 60 degrees to limit the max travel the robot arm has to travel
                
                hexagonPoints = np.transpose(hexagonPoints) # transposes array because it has been being stored in the wrong order
                #print(hexagonPoints)
                hexagonCtr = centerFind(hexagonPoints) #calls centerFind function to evaluate the center of the hexagon
                hexagonCtr[2] = point[2] #stores z value into the z value of the hexagon center array
                
                HexagonCtrs = np.vstack([HexagonCtrs,hexagonCtr]) #adds the hexagon center to the hexagon centers array
                hexcount = 0 #restarts hex count
                hexagonPoints = np.empty((3,5)) #clears hexagon points array 
                HexagonDia = np.vstack([HexagonDia,round(lastmag*np.sqrt(3),3)]) #evaluates heagon diameter with the length of a side
            #print(angle)
            #print(fullHexCount)
        else:
            hmmNewline=np.append(hmmNewline,0) #debug

        #print(hmmNewline)
        hmmArray=np.vstack([hmmArray,hmmNewline]) #debug
        hmmNewline = [] #wipes debuh new line

        #print(delta)
        lastdelta = delta#swaps delta to last delta
        lastmag = curmag #swaps curmag to last mag
        
    
    lastPt = point #swaps point to last point   
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
                #print(hexagonfragNoZ)
                #print(np.shape([hexagonfrag]))
                #print(np.linalg.norm(hexagonfragNoZ-lastFragNoZ))
                #print(np.linalg.norm(fullNuts[i,:]-hexagonfrag))
                if np.linalg.norm(hexagonfragNoZ-CurrentEntryNoZ) <.1:
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
             
fullNuts = fullNuts[2:,:]
np.round(fullNuts,4)
print(fullNuts)
print(np.shape(fullNuts))
        
