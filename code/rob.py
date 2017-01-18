
from math import sqrt, log
import numpy
import gzip
import sys


def read_fasta(fname, whole_id=True):
    """
    Read a fasta file and return a hash. 
    
    If wholeId is set to false only the first part of the ID 
    (upto the first white space) is returned
    """
    try:
        if fname.endswith('.gz'):
            f = gzip.open(fname, 'rb')
        else:
            f = open(fname, 'r')
    except:
        sys.exit("Unable to open file " + fname)

    seqs = {}
    seq = ""
    seqid = ""
    for line in f:
        line = line.rstrip('\r\n')
        if line.startswith(">"):
            if seqid != "":
                seqs[seqid] = seq
                seq = ""
            seqid = line.replace(">", "", 1)
            if not whole_id and seqid.count(" ") > 0:
                seqids = seqid.split(" ")
                seqid = seqids[0]
        else:
            seq += line

    seqs[seqid] = seq
    return seqs


def readFasta(file, whole_id=True):
    """
    Read a fasta file and return a hash.

    If wholeId is set to false only the first part of the ID (upto the first white space) is returned
    """
    return read_fasta(file, whole_id)

def stream_fastq(fqfile):
    """Read a fastq file and provide an iterable of the sequence ID, the
    full header, the sequence, and the quaity scores.

    Note that the sequence ID is the header up until the first space,
    while the header is the whole header.
    """
    
    if fqfile.endswith('.gz'):
        qin = gzip.open(fqfile, 'rb')
    else:
        qin = open(fqfile, 'r')

    while True:
        header = qin.readline()
        if not header:
            break
        header = header.strip()
        seqidparts = header.split(' ')
        seqid = seqidparts[0]
        seq = qin.readline()
        seq = seq.strip()
        qualheader = qin.readline()
        qualscores = qin.readline()
        qualscores = qualscores.strip()
        header = header.replace('@', '', 1)
        yield seqid, header, seq, qualscores


def stdev(arr):
    n = len(arr)
    sm = sum(arr)
    nmpyarr = numpy.array(arr)
    sum_of_squares = numpy.sum(nmpyarr * nmpyarr)
    result = 0
    try:
        result = sqrt(sum_of_squares / n - (sm / n) ** 2)
    except:
        sys.stderr.write("can't create the stdev. Only have " + str(n) + " elements in the array")

    return result


def mean(arr):
    n = len(arr)
    sm = sum(arr)
    return float(sm)/float(n)


def median(arr):
    arr.sort()
    n = len(arr)
    mid = n/2
    return arr[mid]


def rc(dna):
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    rcseq = dna.translate(complements)[::-1]
    return rcseq


def shannon(dna, word):
    count = {}
    for i in range(len(dna)-word):
        substr = dna[i:i+word]
        count[substr] = count.get(substr, 0) + 1

    hp = 0
    p = len(dna) - word
    for s in count:
        n = 1.0 * count[s]
        hp += float((n/p) * log(n/p))
    return 0-hp


def ed_neighborsof(i,j, scorematrix):
    neighbors = []

    if i - 1 >= 0 and j - 1 >= 0:
        neighbors.append(scorematrix[i-1][j-1])
        neighbors.append(scorematrix[i-1][j])
        neighbors.append(scorematrix[i][j-1])
    elif i - 1 >= 0:
        neighbors.append(scorematrix[i-1][j])
        #neighbors.append(j+1)
    elif j - 1 >= 0:
        neighbors.append(scorematrix[i][j-1])
        #neighbors.append(i+1)
    else:
        neighbors.append(0)

    return neighbors

def ed_penalty(acid1, acid2):
    return 1


def edit_distance(seq1, seq2):

    scorematrix = [[0 for s in seq1] for t in seq2]

    for i in range(0, len(seq2)):
        for j in range(0, len(seq1)):
            if seq2[i] == seq1[j]:
                if i > 0 and j > 0:
                    scorematrix[i][j] = scorematrix[i-1][j-1]
                else:
                    scorematrix[i][j]= min(ed_neighborsof(i,j,scorematrix))
            else:
                scorematrix[i][j] = min(ed_neighborsof(i,j,scorematrix)) + ed_penalty(seq2[i], seq1[j])

    return scorematrix[len(seq2) - 1][len(seq1) - 1]





import string


def ascii_clean(s):
    """Remove non-ascii characters from a string"""
    return filter(lambda x: x in string.printable, s)


import math


def latlon2distance(lat1, long1, lat2, long2, miles=False):
    """Convert two coordinates to distance. 

    This is an approximation since the earth is not spherical, but accuracy is <100m, especially for close points
    
    This code was taken from http://www.johndcook.com/python_longitude_latitude.html

    Latitude is measured in degrees north of the equator; southern locations have negative latitude. 
    Similarly, longitude is measured in degrees east of the Prime Meridian. A location 10deg west of 
    the Prime Meridian, for example, could be expressed as either 350deg  east or as -10deg east.

    Arguments: lat1, long1; lat2, long2; miles is a boolean. If you want miles set it to true. Else set it to false

    """


    # Convert latitude and longitude to 
    # spherical coordinates in radians.
    degrees_to_radians = math.pi/180.0

    # phi = 90 - latitude
    phi1 = (90.0 - lat1)*degrees_to_radians
    phi2 = (90.0 - lat2)*degrees_to_radians

    # theta = longitude
    theta1 = long1*degrees_to_radians
    theta2 = long2*degrees_to_radians

    # Compute spherical distance from spherical coordinates.

    # For two locations in spherical coordinates 
    # (1, theta, phi) and (1, theta, phi)
    # cosine( arc length ) = 
    #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length

    cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) + math.cos(phi1)*math.cos(phi2))
    arc = math.acos( cos )

    # Remember to multiply arc by the radius of the earth 
    # in your favorite set of units to get length.
    #
    # To convert to miles multiple arc by 3960
    # To convert to kilometers multiply arc by 6373

    if miles:
        arc *= 3960
    else:
        arc *= 6373

    return arc


# run this as a script
if __name__ == "__main__":
    import sys
