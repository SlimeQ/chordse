import numpy
from scipy.io import wavfile

if __name__ == '__main__':
	sampleRate, data = wavfile.read("im_being_watched_by_the_cia.wav")
	numpy.savetxt("im_being_watched_by_the_cia.csv", data[:,0], delimiter=",", fmt="%10.5f")