'''
List complete bacteria or phage IDs. Set the bool to true to get bacteria
'''
import sys
from phage import Phage

phage=Phage()

try:
    bacteria = sys.argv[1]
except:
    bacteria = False


if bacteria:
    d = phage.completeBacteriaIDs()
else:
    d = phage.phageIDs()

for p in d:
    print(p)
