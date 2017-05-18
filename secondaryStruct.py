class SecondaryStructFASTA(object):
	"""Class used to store a FASTA object that has
	secondary protein structure information.
	
	Attributes:
		header - Contains the Header of the FASTA file
		sequence - Contains the sequence in the FASTA file
		secondaryStruct - Contains the secondary protein structure information
	
	TODO: This code doesn't meet the FASTA specs:
		see <https://en.wikipedia.org/wiki/FASTA_format> for more info
	"""
	def __init__(self, header, sequence, secondaryStruct):
		"""Constructor, fileName is a string showing the fileName to load """
		self.header = header
		self.sequence = sequence
		self.secondaryStruct = secondaryStruct
		
def outputWithoutSecondaryStruct(SecondaryStructFASTAs, fileName):
	""" Writes the SecondaryStructFASTAs without the secondaryStruct data.

		This is so for testing online algorithms.
		
		Args:
			SecondaryStructFASTAs: a list of SecondaryStructFASTA objects.
			fileName: the file to save the data
	"""
	with open(fileName, "w") as f:
		for SecondaryStructFASTA in SecondaryStructFASTAs:
			f.write(">" + SecondaryStructFASTA.header)
			f.write(SecondaryStructFASTA.sequence)

def loadSecondaryStructFASTAs(fileName, limit = 2000, seqLengthLimit = 2000):
	""" Loads a list of FASTA objects with secondary structure information.
	
	The file should have the following format:
	>header
	SEQUENCE
	SECONDARY_STRUCTUE_INFORMATION
	>header2
	... etc
	
	Args:
		fileName: the file-name of the list of FASTA objects.
		limit: The max amount of FASTA objects to load.
		seqLengthLimit: the max length of a sequence (some algorithms have limits)
		
	Returns a list of SecondaryStructFASTA objects.
	"""
	secondaryStructFASTAs = []
	
	with open(fileName, "r") as f:
		numberFound = 0
		while True:
			line = f.readline()
			
			if not line: #finished reading
				return secondaryStructFASTAs
			if line[0] != ">": #not the beginning of a FASTA
				print("Didn't find a header, file may be incorrect.")
				continue
			
			header = line[1:] #remove leading >
			sequence = f.readline()
			secondaryStruct = f.readline()
			
			if(len(sequence) >= seqLengthLimit):
				continue #too long so can't use
				
			secondaryStructFASTA = SecondaryStructFASTA(header, sequence,
			                                            secondaryStruct)
			secondaryStructFASTAs.append(secondaryStructFASTA)
			numberFound+=1
			
			if numberFound >= limit:
				return secondaryStructFASTAs


