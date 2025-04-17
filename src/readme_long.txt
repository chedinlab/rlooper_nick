#Biophysical Parameter

--N [int]
Followed by an integer specifies the size of the superhelical domain in base-pairs. This is the number of bases over which DNA is not in its relaxed state. Alternatively, --N auto can be used to set the size of the superhelical domain to the size of the provided sequence. The default if unspecified is 1500bp (Kouzine and Levins, Myc C FUSE work).


--sigma [float]
Is the superhelical density over a region of size N. The default value is -0.07, or 7% negative supercoiling, which is the physiological level of supercoiling in many Prokaryotes, and also a rough approximation of the amount of supercoiling being driven by transcription in Eukaryotes.


--a [float]
Followed by an floating point value specifies the energetic disfavorability in kcal/mol associated with the two junctions at each end of the R-loop structure. This quantity is assumed to be length independent. The default value if unspecified is 10 kcal/mol (5 per junction), consistent with the experimentally measured value for other types of superhelically driven DNA structural transition. (Levins, Benham)


--minlength [int]
Followed by an integer specifies the size of R-loops in base-pairs below which should be excluded from the biophysical ensemble. The default if unspecified is 2bp, which is the full ensemble.


#Sequence Handling Overrides

--reverse [bool]
Reverses the transcribed orientation of the provided input sequence. Required for sequence provided in the 3' to 5' orientation.


--complement [bool]
Complements the transcribed orientation of the provided input sequence.
