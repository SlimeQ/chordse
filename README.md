# chordse

This application takes input in the form of an audio file and outputs vectors representing the presence of each note in the musical scale. The resulting vectors can then be used for chord recognition tasks.

chordse.py is the serial version of this program.
chordsepool.py is a variation which runs in parallel across a number of threads.
chromagram.c is a C module used by the python scripts in order to improve speed.
