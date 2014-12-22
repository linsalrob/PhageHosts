import os
import sys


try:
    phageIdF=sys.argv[1]
    bactIdF=sys.argv[2]
    predF=sys.argv[3]
except:
    sys.exit(sys.argv[0] + " <phage id file> <bacteria id file> <predictions file>\nThe id files should have [id, genus species]. The prediction file should have [phage id, score, predicted host(s)] with multiple equal predictions separated by tabs\n")


phage={}
with open(phageIdF, 'r') as fin:
    for line in fin:
        line=line.strip()
        p=line.split("\t")
        phage[p[0]]=p[1]

bact={}
with open(bactIdF, 'r') as fin:
    for line in fin:
        line=line.strip()
        p=line.split("\t")
        bact[p[0]]=p[1]


print "\t".join(["Phage", "Real host", "Predicted host"])
same=0
diff=0
with open(predF, 'r') as pred:
    for line in pred:
        line=line.strip()
        p=line.split("\t")
        p=[x.strip() for x in p]
        if p[0] not in phage:
            sys.stderr.write("No phage was found for |" + str(p[0]) + "|\n")
            continue
        real = phage[p[0]]

        host={}
        for h in p[2:len(p)]:
            if h in bact:
                if bact[h] in host:
                    host[bact[h]] += 1
                else:
                    host[bact[h]] = 1
            else:
                sys.stderr.write("No host was found for " + str(h) + "\n")

        prediction = None
        if len(host) == 1:
            prediction = host.keys()[0]
        else:
            most = 0
            for h in sorted(host, key=host.get, reverse=True):
                if host[h] > most:
                    most=host[h]
                    prediction = h
                elif host[h] == most and h == real:
                    prediction = h

        if prediction == real:
            same += 1
        else:
            diff += 1

        print "\t".join([p[0], real, prediction])
            
sys.stderr.write("Overall: correct predictions: " + str(same) + " incorrect predictions " + str(diff) + "\n")
