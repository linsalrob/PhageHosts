import os
import sys

'''

The ROC table looks like this:

    From            To              Distance        speci   genus   family  order   class   phylum  superkingdom
    NC_024375       NC_005042       0.964027348451  0       0       0       0       0       0       1
    NC_024375       NC_005043       0.956167698873  0       0       0       0       0       0       1
    NC_024375       NC_005823       0.935176712175  0       0       0       0       0       0       1

We need to find the lowest score in column 2 and compare the counts

'''

try:
    f = sys.argv[1]
except:
    sys.exit(sys.argv[0] + " <ROC output file>")

data = {}
with open(f, 'r') as fin:
    for l in fin:
        p=l.strip().split("\t")
        p1 = p.pop(0)
        p2 = p.pop(0)

        p[0] = float(p[0])
        for i in range(1, len(p)):
            p[i] = int(p[i])
        

        if p1 not in data:
            data[p1] = [p2] + p
            continue

        if p[0] > data[p1][1]:
            #sys.stderr.write("Had " + data[p1][0] + " (" + str(data[p1][1])+ ") and replaced it with " + p2 + " (" + str(p[0]) + ")\n") 
            data[p1] = [p2] + p
        elif p[0] == data[p1][1]:
            oldsum = sum(data[p1][2:])
            newsum = sum(p[1:])
            #sys.stderr.write("EQUAL:\n" + p1 + "\t" + "\t".join(map(str, data[p1])) + "\n" + p1 + "\t" + "\t".join(map(str, [p2] + p)) + "\n")
            if newsum > oldsum:
                data[p1] = [p2] + p
        
        #sys.stderr.write(str(p[0]) + " ne " + str(data[p1][1]) + "\n")


count = [0, 0, 0, 0, 0, 0, 0]
total = 0


sys.stderr.write("There are " + str(len(data.keys())) + " elements in data\n")

for p in data:
    #sys.stderr.write("BEst hit: " + p + "\t" + "\t".join(map(str, data[p])) + "\n")
    total += 1
    for i in range(len(count)):
        count[i] += data[p][i+2]

print("Taxonomic level\tCorrect assignments\tIncorrect assignments\tTotal Assignments\tPercent Correct")
taxonomy = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
for i in range(len(taxonomy)):
    fr = "%0.2f" % (1.0 * count[i]/total*100)
    print(taxonomy[i] + "\t" + str(count[i]) + "\t" + str(total - count[i]) + "\t" + str(total) + "\t" + fr)



           





