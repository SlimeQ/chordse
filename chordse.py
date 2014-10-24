from __future__ import division
from scipy.io import wavfile
import numpy as np
import numpy.ctypeslib as npct
import math
import ctypes
from datetime import datetime

# B = bins
B = 12
fmin = 96
fmax = 5250
M = int(math.log(fmax/fmin, 2))
K = int(math.ceil(B * M))
Q = 1 / (2**(1/B) - 1)
numChunks = 500
i = 0

def chromagram(Xcq):
    global B
    print "chromagramming"
    # CH = [sum([math.fabs(Xcq[b + m * B]) for m in xrange(M)]) for b in xrange(B)]
    

    CH = []
    for b in xrange(B):
        summation = 0
        for m in xrange(M):
            summation += math.fabs(Xcq[b+m*B])
        CH.append(summation)
    CH = np.subtract(CH, np.min(CH))
    CH = np.divide(CH, np.max(CH))
    print CH
    # import code
    # code.interact(local=locals())
    return CH

def constantQ(fs, x):
    print "constant Q-ing"
    # global i
    # fs = sampling period
    # x = data

    # averaging the two channels for some reason
    x = [y[0] for y in x] # np.mean(x, axis=1)

    X = []
    for k in xrange(K):
        fk = 2**(k/B) * fmin
        N = min(int(math.ceil((Q * fs) /fk)), len(x))
        W = np.hamming(N)

        summation = 0
        for n in xrange(N):
        #     print n, ":", N
        #     print "len(x) =", len(x)
        #     print "len(W) =", len(W)
            summation += x[n] * W[n] * np.exp(-2 * np.pi * 1j * Q * n / N)

        X.append(summation / N)

    # print X
    # print i
    # i += 1
    return X

    



if __name__ == '__main__':

    # read file
    sampleRate, data = wavfile.read("im_being_watched_by_the_cia.wav")

    # chunk datas into numChunks chunks of N length
    chunks = np.array_split(data, numChunks)
    N = len(chunks[0])

    # link to c module
    lib = ctypes.CDLL("./chromagram.so")
    chromaFunc = lib.chromagram
    multiChroma = lib.multichroma


    out = np.zeros((B * numChunks), dtype=np.float64)
    # for x in data[:,0]:
    #     print x,
    multiChroma(npct.as_ctypes(np.ascontiguousarray(data[:,0], dtype=np.int16)), len(data), numChunks, N, fmin, fmax, M, B, K, ctypes.c_float(Q), sampleRate, npct.as_ctypes(np.ascontiguousarray(out, dtype=np.float64)))


    t1 = datetime.now()
    chroma = []
    for chunk in chunks:
        out = np.zeros((B), dtype=np.float64)
        # print chunk[:,0]
        # print max(chunk[:,0])

        chromaFunc(npct.as_ctypes(np.ascontiguousarray(chunk[:,0], dtype=np.int16)), N, fmin, fmax, M, B, K, ctypes.c_float(Q), sampleRate, npct.as_ctypes(np.ascontiguousarray(out, dtype=np.float64)))
        chroma.append(out)
    dt = datetime.now() - t1
    ms = (dt.days * 24 * 60 * 60 + dt.seconds) * 1000 + dt.microseconds / 1000.0
    print "Total time =", ms, "ms"
    # print chroma

  

    # for chunk in chroma:
    #     print chunk
    # print N
