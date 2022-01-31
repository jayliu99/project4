# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
	""" Class for NeedlemanWunsch Alignment

	Parameters:
		sub_matrix_file: str
			Path/filename of substitution matrix
		gap_open: float
			Gap opening penalty
		gap_extend: float
			Gap extension penalty

	Attributes:
		seqA_align: str
			seqA alignment
		seqB_align: str
			seqB alignment
		alignment_score: float
			Score of alignment from algorithm
		gap_open: float
			Gap opening penalty
		gap_extend: float
			Gap extension penalty
	"""
	def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
		# Init alignment and gap matrices
		self._align_matrix = None
		self._gapA_matrix = None
		self._gapB_matrix = None

		# Init matrices for backtrace procedure
		self._back = None
		# self._back_A = None
		# self._back_B = None

		# Init alignment_score
		self.alignment_score = 0

		# Init empty alignment attributes
		self.seqA_align = ""
		self.seqB_align = ""

		# Init empty sequences
		self._seqA = ""
		self._seqB = ""

		# Setting gap open and gap extension penalties
		self.gap_open = gap_open
		assert gap_open < 0, "Gap opening penalty must be negative."
		self.gap_extend = gap_extend
		assert gap_extend < 0, "Gap extension penalty must be negative."

		# Generating substitution matrix
		self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

	def _read_sub_matrix(self, sub_matrix_file):
		"""
		DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

		This function reads in a scoring matrix from any matrix like file.
		Where there is a line of the residues followed by substitution matrix.
		This file also saves the alphabet list attribute.

		Parameters:
			sub_matrix_file: str
				Name (and associated path if not in current working directory)
				of the matrix file that contains the scoring matrix.

		Returns:
			dict_sub: dict
				Substitution matrix dictionary with tuple of the two residues as
				the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
		"""
		with open(sub_matrix_file, 'r') as f:
			dict_sub = {}  # Dictionary for storing scores from sub matrix
			residue_list = []  # For storing residue list
			start = False  # trigger for reading in score values
			res_2 = 0  # used for generating substitution matrix
			# reading file line by line
			for line_num, line in enumerate(f):
				# Reading in residue list
				if '#' not in line.strip() and start is False:
					residue_list = [k for k in line.strip().upper().split(' ') if k != '']
					start = True
				# Generating substitution scoring dictionary
				elif start is True and res_2 < len(residue_list):
					line = [k for k in line.strip().split(' ') if k != '']
					# reading in line by line to create substitution dictionary
					assert len(residue_list) == len(line), "Score line should be same length as residue list"
					for res_1 in range(len(line)):
						dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
					res_2 += 1
				elif start is True and res_2 == len(residue_list):
					break
		return dict_sub

	def _fill_in_mat(self, matrix_to_fill, r, c):
		"""
		This function handles the details of filling in the four matrices
		(3 scoring, 1 backtrace) used to implement the Needleman-Wunsch Algorithm.

		Parameters:
			matrix_to_fill: char
				Identifies the matrix to be filled in 
				(should be 'M', 'X', 'Y', or 'B')
			r: int
				Identifies the row index of the matrix to be filled in 
			c: int
				Identifies the column index of the matrix to be filled in 
		"""

		assert matrix_to_fill in ['M', 'X', 'Y', 'B'], "Matrix should be M, X, Y, or B"

		# Fill in M matrix
		if matrix_to_fill == 'M':
			options = [self._align_matrix[r-1, c-1],
					   self._gapA_matrix[r-1, c-1],
					   self._gapB_matrix[r-1, c-1]]

			optimal_move = max(options)
			match_score = self.sub_dict[(self._seqA[c-1], self._seqB[r-1])]
			self._align_matrix[r, c] = match_score + optimal_move

		# Fill in X matrix
		elif matrix_to_fill == 'X':
			options = [self.gap_open + self.gap_extend + self._align_matrix[r-1, c],
					   self.gap_extend + self._gapA_matrix[r-1, c],
					   self.gap_open + self.gap_extend + self._gapB_matrix[r-1, c]]

			optimal_move = max(options)
			self._gapA_matrix[r, c] = optimal_move

		# Fill in Y matrix
		elif matrix_to_fill == 'Y': 
			options = [self.gap_open + self.gap_extend + self._align_matrix[r, c-1],
					   self.gap_open + self.gap_extend + self._gapA_matrix[r, c-1],
					   self.gap_extend + self._gapB_matrix[r, c-1]]

			optimal_move = max(options)
			self._gapB_matrix[r, c] = optimal_move

		else: # Fill in backtracking matrix
			options = [self._align_matrix[r,c], # Score if match 
						self._gapA_matrix[r,c],	# Score if insert a gap in seqA
						self._gapB_matrix[r,c]] # Score if insert a gap in seqB
			optimal_move = max(options)
			self._back[r,c] = options.index(optimal_move)


	def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
		"""
		This function performs the Needleman-Wunsch Algorithm for 
		global sequence alignment of seqA and seqB and returns the optimal
		alignment along with the optimal alignment score.

		Parameters:
			seqA: str
				String of characters representing the first sequence to be aligned
				(Also called: Sequence X)
			seqB: str
				String of characters representing the first sequence to be aligned
				(Also called: Sequence Y)

		Returns:
			self._backtrace(): tuple
				A tuple with the following format: 
				(alignment score, seqA alignment, seqB alignment) 
				e.g. (17, "MAVHQLIRRP", "M---QLIRHP")
		"""
		# Initialize matrix private attributes for use in alignment
		# create matrices for alignment scores and gaps
		self._align_matrix = np.ones((len(seqB) + 1, len(seqA) + 1)) * -np.inf
		self._gapA_matrix = np.ones((len(seqB) + 1, len(seqA) + 1)) * -np.inf
		self._gapB_matrix = np.ones((len(seqB) + 1, len(seqA) + 1)) * -np.inf

		# create matrix used in backtrace procedure
		self._back = np.ones((len(seqB) + 1, len(seqA) + 1)) * -np.inf

		# Resetting alignment in case method is called more than once
		self.seqA_align = ""
		self.seqB_align = ""

		# Resetting alignment score in case method is called more than once
		self.alignment_score = 0

		# Initializing sequences for use in backtrace method
		self._seqA = seqA
		self._seqB = seqB

		# Fill out base case for M matrix
		self._align_matrix[0, 0] = 0 

		# Fill out base case for X matrix
		self._gapA_matrix[0, 0] = self.gap_open
		for r in range(1, self._gapA_matrix.shape[0]):
			self._gapA_matrix[r, 0] = self.gap_open + self.gap_extend * r

		# Fill out base case for Y matrix
		self._gapB_matrix[0, 0] = self.gap_open
		for c in range(1, self._gapB_matrix.shape[1]):
			self._gapB_matrix[0, c] = self.gap_open + self.gap_extend * c

		# Fill out base case for backtracing matrix
		for r in range(1, self._back.shape[0]):
			self._back[r, 0] = '1' # If seqA has run out of letters to match, insert gaps until seqB is done
		for c in range(1, self._back.shape[1]):
			self._back[0, c] = '2' # If seqB has run out of letters to match, insert gaps until seqA is done

		# Fill in rest of matrices from left->right, top->bottom
		for r in range(1, self._align_matrix.shape[0]):
			for c in range(1, self._align_matrix.shape[1]):
				self._fill_in_mat('M', r, c)
				self._fill_in_mat('X', r, c)
				self._fill_in_mat('Y', r, c)
				# Fill in backtracking matrix
				self._fill_in_mat('B', r, c)

		return self._backtrace()

	def _do_alignment(self, align_A, align_B, indexA, indexB, r, c) -> Tuple[int, int, int, int]:
		"""

		This function handles the details of filling in the final global alignment
		of seqA and seqB. It is a helper function, meant to be called by a 
		function that executes a backtracing heuristic compatible with Needleman-Wunsch.

		In addition to filling out the seqA and seqB alignments, _do_alignment returns a 
		tuple containing updated indices (to keep track of which letters in seqA and seqB 
		have yet to be aligned) as well as updated row and column values (to keep track of
		which cell in the backtracking matrix contains the next best step).

		Parameters:
			align_A: list
				Maintains seqA's optimal sequence alignment
			align_B: list
				Maintains seqB's optimal sequence alignment
			indexA: int
				Marks the current character to be aligned in seqA 
			indexB: int
				Marks the current character to be aligned in seqB
			r: int
				Identifies the row index of the cell in the
				backtracking matrix currently being considered.
			c: int
				Identifies the column index of the cell in the
				backtracking matrix currently being considered.

		Returns:
			Tuple[int, int, int, int]
				A tuple with the following format: 
				(updated index for seqA, updated index for seq B, updated row, updated col) 
		"""

		if self._back[r,c] == 0: # Match
			align_A.append(self._seqA[indexA])
			align_B.append(self._seqB[indexB])
			return (indexA-1, indexB-1, r-1, c-1)

		elif self._back[r,c] == 1: # Gap in A
			align_A.append("-")
			align_B.append(self._seqB[indexB])
			return (indexA, indexB-1, r-1, c)

		else: # Gap in B
			align_A.append(self._seqA[indexA])
			align_B.append("-")
			return (indexA-1, indexB, r, c-1)


	def _backtrace(self) -> Tuple[float, str, str]:
		"""
		This function performs a traceback procedure based on the
		heuristic used to implement the Needleman-Wunsch Algorithm. 
		It returns the optimal alignment along with the optimal alignment score.

		Returns:
			result: tuple
				A tuple with the following format: 
				(alignment score, seqA alignment, seqB alignment) 
				e.g. (17, "MAVHQLIRRP", "M---QLIRHP")

		"""
		# Implement this method based upon the heuristic chosen in the align method above.
		
		align_A = []
		align_B = []

		indexA = len(self._seqA) - 1
		indexB = len(self._seqB) - 1

		r = self._back.shape[0] - 1
		c = self._back.shape[1] - 1

		# Perform backtracking
		while (r,c) != (0,0):
			indexA, indexB, r, c = self._do_alignment(align_A, align_B, indexA, indexB, r, c)

		# Calcuate alignment score
		score = 0
		gap_start = True
		align_A.reverse()
		align_B.reverse()
		for i in range(len(align_A)):
			if align_A[i] == "-" or align_B[i] == "-":
				if gap_start:
					score += self.gap_open + self.gap_extend
					gap_start = False
				else: 
					score += self.gap_extend
			else: 
				score += self.sub_dict[align_A[i], align_B[i]]
				gap_start = True

		# Output final aligned sequences and alignment score
		result = (score, ''.join(align_A), ''.join(align_B))
		return result



def read_fasta(fasta_file: str) -> Tuple[str, str]:
	"""
	DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

	This function reads in a FASTA file and returns the associated
	string of characters (residues or nucleotides) and the header.
	This function assumes a single protein or nucleotide sequence
	per fasta file and will only read in the first sequence in the
	file if multiple are provided.

	Parameters:
		fasta_file: str
			name (and associated path if not in current working directory)
			of the Fasta file.

	Returns:
		seq: str
			String of characters from FASTA file
		header: str
			Fasta header
	"""
	assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
	with open(fasta_file) as f:
		seq = ""  # initializing sequence
		first_header = True
		for line in f:
			is_header = line.strip().startswith(">")
			# Reading in the first header
			if is_header and first_header:
				header = line.strip()  # reading in fasta header
				first_header = False
			# Reading in the sequence line by line
			elif not is_header:
				seq += line.strip().upper()  # generating full sequence
			# Breaking if more than one header is provided in the fasta file
			elif is_header and not first_header:
				break
	return seq, header
