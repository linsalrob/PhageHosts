'''
Mutual information 

Suggested by Peter Salamon:
    You need a threshold of present in the metagenome to turn your
    current variables into 0,1 variables that signal presence of the
    entity in the metagenome. As a start, you could threshold on mean
    value for that entity (phage or bacterium). Then you need the
    fraction of metagenomes in which the phage and the bacterium 
    co-occur, the fraction with one not the other and the fraction with
    the other and not the one and the fraction with neither. This gives
    the joint distribution of the four cases. You also need the marginal
    distributions of fraction in which the phage occurs and the fraction
    in which the bacterium occurs. Then simply calculate the Kullback
    entropy of the joint distribution relative to the product
    distribution.

In Peter's notes (see email Dec 16th):
    P1 = Z/3205 = phagebacteriasum / nmg
    P2 = (X - Z) / 3205 = (bacteriasum - phagebacteriasum) / nmg
    P3 = (Y - Z) / 3205 = (phagesum - phagebacteriasum) / nmg 
    P4 = 1 - P1 - P2 - P3

    Q1 = XY / 3205^2  = (phagesum * bacteriasum) / nmg**2
    Q2 = (3205 - X)Y / 3205^2 = ((nmg - bacteriasum) * phagesum ) / nmg**2
    Q3 = X(3205-Y) / 3205^2 = (bacteriasum * (nmg-phagesum)) / nmg**2
    Q4 = 1 - Q1 - Q2 - Q3


'''

import sys,os
import math
import numpy

# we need to get rid of the nans
def notnan(x): return not math.isnan(x)
# and this to get rid of the zeros
def notzero(x): return x > 0

# a function for the stdev
def stdev(arr):
    n = len(arr)
    sm = sum(arr)
    nmpyarr = numpy.array(arr)
    sum_of_squares = numpy.sum(nmpyarr * nmpyarr)
    result=0
    try:
        result = math.sqrt(sum_of_squares / n - (sm / n) ** 2)
    except:
        sys.stderr.write("There was an error trying to calculate the stdev. sum = " + str(sm) + "; n = " + str(n) + "; ss = " + str(sum_of_squares) + "\n")

    return result



try:
    phageF = sys.argv[1]
    bacteriaF = sys.argv[2]
    stdevs = sys.argv[3]
except:
    sys.exit(sys.argv[0] + " <phage metagenome presence>  <bacterial metagenome presence> <number of stdevs below mean to include as positives. Use 0 to use the mean as the cutoff. Use -1 to ignore all 0 values in the data (and use the mean as a cutoff)>; use -2 to use any value to indicate present")

stdevs=float(stdevs)
ignoreZeros = False
any_value = False
if stdevs == -1:
    ignoreZeros = True
    stdevs = 0

if stdevs == -2:
    any_value = True

phage={}
with open(phageF, 'r') as pin:
    l=pin.readline()
    headers=l.strip().split("\t")
    headers.pop(0)
    headers.pop()
    nmg = len(headers)
    for l in pin:
        p = l.strip().split("\t")
        phageId = p.pop(0)
        phage[phageId]={}
        # remove the taxonomy
        p.pop()
        p=map(float, p)
        data = filter(notnan, p)
        if ignoreZeros:
            data=filter(notzero, data)
        if len(data) == 0:
            for i in xrange(len(p)):
                phage[phageId][headers[i]]=0
            continue
        mean = sum(data)/len(data)
        s = stdev(data) * stdevs
        threshold = mean - s
        # sys.stderr.write("Phage: " + phageId + " threshold " + str(threshold) + " mean=" + str(mean) + " s="  + str(s) + "\n")
        if any_value:
            threshold = 0
        for i in xrange(len(p)):
            if math.isnan(p[i]) or p[i] <= threshold:
                phage[phageId][headers[i]]=0
            else:
                phage[phageId][headers[i]]=1



bacteria={}
with open(bacteriaF, 'r') as pin:
    l=pin.readline()
    headers=l.strip().split("\t")
    headers.pop(0)
    headers.pop()
    assert len(headers) == nmg
    for l in pin:
        p = l.strip().split("\t")
        bacteriaId = p.pop(0)
        bacteria[bacteriaId]={}
        # remove the taxonomy
        p.pop()
        p=map(float, p)
        data = filter(notnan, p)
        if ignoreZeros:
            data=filter(notzero, data)
        if len(data) == 0:
            for i in xrange(len(p)):
                bacteria[bacteriaId][headers[i]]=0
            continue

        mean = sum(data)/len(data)
        s = stdev(data) * stdevs
        threshold = mean - s
        # sys.stderr.write("Bacteria: " + bacteriaId + " threshold " + str(threshold) + " mean=" + str(mean) + " s=" + str(s) + "\n")
        if any_value:
            threshold = 0
        for i in xrange(len(p)):
            if math.isnan(p[i]) or p[i] <= threshold:
                bacteria[bacteriaId][headers[i]]=0
            else:
                bacteria[bacteriaId][headers[i]]=1



allbacteria = bacteria.keys()
print("Phage\tBacteria\tMutual Information")
for p in phage:
    phagesum=0
    for h in phage[p]:
        phagesum += phage[p][h]
    if phagesum == nmg or phagesum == 0:
        for a in allbacteria:
            print(p + "\t" + a + "\t0")
        # these samples are not informative since everything is everywhere or nowhere
        continue
    
    for a in allbacteria:
        bacteriasum=0
        phagebacteriasum=0
        for h in phage[p]:
            bacteriasum+=bacteria[a][h]
            phagebacteriasum += (bacteria[a][h] * phage[p][h])
        # check that the phage and bacteria sum is less than or equal to each individual sum
        assert phagebacteriasum <= phagesum
        assert phagebacteriasum <= bacteriasum

        if bacteriasum == nmg or bacteriasum == 0:
            # these samples are not informative since everything is everywhere or nowhere
            print(p + "\t" + a + "\t0") 
            continue
        
        if phagebacteriasum == 0:
            # these samples never appear together
            print(p + "\t" + a + "\t0") 
            continue

        p1 = 1.0 * phagebacteriasum / nmg
        p2 = 1.0 * (bacteriasum - phagebacteriasum) / nmg
        p3 = 1.0 * (phagesum - phagebacteriasum) / nmg 
        p4 = 1.0 - p1 - p2 - p3

        q1 = 1.0 * (phagesum * bacteriasum) / nmg**2
        q2 = 1.0 * ((nmg - bacteriasum) * phagesum ) / nmg**2
        q3 = 1.0 * (bacteriasum * (nmg-phagesum)) / nmg**2
        q4 = 1.0 - q1 - q2 - q3


        # now we just check for any zero values and ignore those
        if p1 == 0:
            # should not happen!
            v1 = 0
        else:
            v1 = (1.0 * p1 * math.log(p1/q1))

        if p2 == 0:
            # bacteriasum == phagebacteriasum
            v2 = 0
        else:
            v2 = (1.0 * p2 * math.log(p2/q2))

        if p3 == 0:
            # phagesum == phagebacteriasum
            v3 = 0
        else:
            v3 = (1.0 * p3 * math.log(p3/q3))

        if p4 == 0:
            # should not happen!
            v4 = 0
        else:
            v4 = (1.0 * p4 * math.log(p4/q4))

        if q1 == 0:
            sys.stderr.write("ERROR: Q1 is zero")
            sys.stderr.write(" because nmg= " + str(nmg) + "; bacteriasum= " + str(bacteriasum) + "; phagesum= " + str(phagesum) + "; phagebacteriasum=" + str(phagebacteriasum) + "\n")
            print(p + "\t" + a + "\t0") 
            continue

        if q2 == 0:
            sys.stderr.write("ERROR: Q2 is zero")
            sys.stderr.write(" because nmg= " + str(nmg) + "; bacteriasum= " + str(bacteriasum) + "; phagesum= " + str(phagesum) + "; phagebacteriasum=" + str(phagebacteriasum) + "\n")
            print(p + "\t" + a + "\t0") 
            continue

        if q3 == 0:
            sys.stderr.write("ERROR: Q3 is zero")
            sys.stderr.write(" because nmg= " + str(nmg) + "; bacteriasum= " + str(bacteriasum) + "; phagesum= " + str(phagesum) + "; phagebacteriasum=" + str(phagebacteriasum) + "\n")
            print(p + "\t" + a + "\t0") 
            continue

        if q4 == 0:
            sys.stderr.write("ERROR: Q4 is zero")
            sys.stderr.write(" because nmg= " + str(nmg) + "; bacteriasum= " + str(bacteriasum) + "; phagesum= " + str(phagesum) + "; phagebacteriasum=" + str(phagebacteriasum) + "\n")
            print(p + "\t" + a + "\t0") 
            continue
        #sys.stderr.write("p1="+str(p1)+"; p2="+str(p2)+"; p3="+str(p3)+"; p4="+str(p4)+"; q1="+str(q1)+"; q2="+str(q2)+"; q3="+str(q3)+"; q4="+str(q4)+";\n")
        
        
        #mutualinformation = (p1 * math.log(p1/q1)) + (p2 * math.log(p2/q2)) + (p3 * math.log(p3/q3)) + (p4 * math.log(p4/q4))
        mutualinformation = v1 + v2 + v3 + v4
        print(p + "\t" + a + "\t" + str(mutualinformation)) 






