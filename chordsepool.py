from __future__ import division
from scipy.io import wavfile
import scipy.signal as sig
import numpy as np
import numpy.ctypeslib as npct
import math
import ctypes
from multiprocessing import Pool
from datetime import datetime

# B = bins
B = 36
exp2_1_B = ctypes.c_double(2**(1/B))

fmin = 16.35
fmin_ctype = ctypes.c_double(fmin)
fmax = 7902.13
fmax_ctype=ctypes.c_double(fmax)
M = int(math.ceil(math.log(fmax/fmin, 2)))
K = int(B * M)
Q = 1 / (2**(1/B) - 1)
Q_ctype = ctypes.c_double(Q)
target_N = 11025
i = 0

# link to c module
lib = ctypes.CDLL("./chromagram.so")
chromaFunc = lib.chromagram
chromaFunc.restype = ctypes.POINTER(ctypes.c_double)
    
def chromaSuperWrapper(args):
    x = chromaLesserWrapper(*args)
    Array = B*ctypes.c_double 
    x = Array.from_address(ctypes.addressof(x.contents)) 
    x = npct.as_array(x)
    return x

def chromaLesserWrapper(chunk, N, sampleRate, fs):
    # downsample to 11025 Hz
    rateChange = fs / sampleRate
    sig.resample(chunk, len(chunk) * rateChange)
    return chromaFunc(npct.as_ctypes(np.ascontiguousarray(chunk[:,0], dtype=np.int16)), N, fmin_ctype, fmax_ctype, M, B, exp2_1_B, K, Q_ctype, fs)


if __name__ == '__main__':

    # read file
    sampleRate, data = wavfile.read("Eb2.wav")
    # sampleRate, data = wavfile.read("Eb2.wav")
    fs = 11025
    # target_N = fs

    # print "sample rate =", fs

    t1 = datetime.now()

    # chunk datas into C chunks of ~N length
    chunks = np.array_split(data, len(data) / target_N)
    N = len(chunks[0])

    # spawn process pool
    pool = Pool(processes=16)
    chroma = pool.map(chromaSuperWrapper, [(chunk, N, sampleRate, fs) for chunk in chunks])

    dt = datetime.now() - t1
    ms = (dt.days * 24 * 60 * 60 + dt.seconds) * 1000 + dt.microseconds / 1000.0
    
    for vector in chroma:
        print '(' + ', '.join([str(i) for i in vector]) + ')'
        import code
        code.interact(local=locals())

    print "Total time =", ms, "ms"
    print "------------\n"

        
    exit()  
