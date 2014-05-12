@outputSchema("genesWithBox: {(chr_f:chararray, strand_f:chararray, start_f:int, end_f:int, name_f:chararray, box:int)}")
def foo(chr_f, strand_f, start_f, end_f, name_f, box):
	out = []
	nBoxes = end_f/box - start_f/box 
	while nBoxes > -1:
		out.append(tuple([chr_f, strand_f, start_f ,end_f ,name_f, nBoxes]))
		nBoxes -=1 
	return out
