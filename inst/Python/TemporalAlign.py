import numpy as np
import pandas as pd
import math as math
import re

pd.options.display.max_columns = None

#Initialise the H matrix, or score matrix, that will be used for scoring according to sequence size
def init_Hmat(s1_len,s2_len,g,local_align):

	#For a global alignment, we want to penalise gaps at the beginning and end of each sequence by an incrementing gap penalty
	if local_align == 0:
		H = np.zeros((s2_len+1,s1_len+1), float)

		#Fill the first column with incrementing gap penalties
		for row in range(s2_len+1):
			H[row][0] = -row*g
	
		#Fill the first row with incrementing gap penalties
		for col in range(s1_len+1):
			H[0][col] = -col*g
	
		return H

	#For a local alignment, we do not want to penalise gaps at the beginning and end of each sequence
	elif local_align == 1:
		H = np.zeros((s2_len+1,s1_len+1), float)

		#Fill the first column with incrementing 0's
		for row in range(s2_len+1):
			H[row][0] = 0
	
		#Fill the first row with 0's
		for col in range(s1_len+1):
			H[0][col] = 0
	
		return H

#Initialise the TR matrix that will track the row time penalty
def init_TRmat(s1,s1_len,s2,s2_len,local_align):
	TR = np.zeros((s2_len+1,s1_len+1), float)

	#Fill in a time penalty depending on the time gap stored in the corresponding seq2 vector
	for row in range(1,s2_len+1):
		TR[row][0] = TR[row-1][0] + float(s2[row-1][0])

	return TR

#Initialise the TC matrix that will track the row time penalty
def init_TCmat(s1,s1_len,s2,s2_len,local_align):
	TC = np.zeros((s2_len+1,s1_len+1), float)

	#Fill in a time penalty depending on the time gap stored in the corresponding seq2 vector
	for col in range(1,s1_len+1):
		TC[0][col] = TC[0][col-1] + float(s1[col-1][0])

	return TC

#Initialise the traceback matrix that will track the direction we took in H
def init_traceMat(s1_len,s2_len,local_align):

	#In the case of global alignment, we will need to progress to the top-left square of our trace matrix, and fill the top row
	#and col with the according direction, 2 being left and 1 being up.
	if local_align == 0:
		traceMat = np.zeros((s2_len+1,s1_len+1), float)
		#Fill top row with 2
		traceMat[0] = 2
		#Fill left col with 1
		traceMat[:,0] = 1
		#Fill in DONE square
		traceMat[0][0] = -1
	
	#For a local alignment, we simply wish to end our alignment traceback if we reach the end of a sequence, and thus we put
	#a 0 in each relevant cell.
	elif local_align == 1:
		traceMat = np.zeros((s2_len+1,s1_len+1), float)
		#Fill top row with 0
		traceMat[0] = 0
		#Fill left col with 0
		traceMat[:,0] = 0
		#Fill in DONE square
		traceMat[0][0] = 0

	return traceMat

#Input the drugs x and y and returns their corresponding cell value in the subsitution matrix, s
def score(s,x,y):
	score = s[x][y]

	return score

#Temporal Needleman Wunsch score matrix calculation
#Here we fill our H matrix with values, starting from the top left, according to those equations defined in 10.1109/dsaa.2015.7344785
def TNW_scoreMat(s1,s1_len,s2,s2_len,g,T,H,TR,TC,traceMat,s):
	#initialise time penalty at 0
	tp = 0

	#Initialise final alignment score#
	finalScore = 0

	for i in range(1,s2_len+1):
		for j in range(1,s1_len+1):

			if i == 1 and j == 1:
				tp = 0
			#Calculate correct time penalty according to ref.
			else:
				tpx = float(s2[i-1][0]) + float(TR[i-1][j-1])
				tpy = float(s1[j-1][0]) + float(TC[i-1][j-1])

				tpmax = max(tpx,tpy)

				#Avoid division by 0
				if tpmax == 0:
					tpmax = 1e-8
				
				tp = T*(abs(tpx-tpy)/tpmax)

			#Calculate potential Hmat choices
			score_match = score(s,s1[j-1][1],s2[i-1][1])
			Hup = H[i-1][j-1] + score_match - tp
			Hmid = H[i-1][j] - g
			Hbot = H[i][j-1] - g

			#Select absolute maximum of options
			H[i][j] = max([Hup,Hmid,Hbot])

			finalScore += H[i][j]

			#Update traceMat
			traceMat[i][j] = [Hup,Hmid,Hbot].index(H[i][j])
			traceVal = traceMat[i][j]

			#Update TR and TC matrices according to choice
			if traceVal == 0:
				TR[i][j] = 0
				TC[i][j] = 0
			elif traceVal == 1:
				TR[i][j] = TR[i-1][j] + float(s2[i-1][0])
				TC[i][j] = TC[i-1][j]
			elif traceVal == 2:
				TR[i][j] = TC[i][j-1]
				TC[i][j] = TC[i][j-1] + float(s1[j-1][0])

	return finalScore

#Here we utilise our traceback matrix to fill in all matches and gaps in each sequence to return the aligned sequences
#which achieved the highest scoring path in H
def align_TNW(traceMat, s1, s2, s1_len, s2_len):
	s1_aligned = ""
	s2_aligned = ""

	#Initiate countdown variables for lengths of s1 and s2
	j = s1_len
	i = s2_len

	#Count total pairs of perfect alignment
	totAligned = 0

	while j > 0 and i > 0:
		#Perfect match, diagonal move
		if traceMat[i][j] == 0:
			s1_aligned = s1[j-1][1] + s1_aligned
			s2_aligned = s2[i-1][1] + s2_aligned
			i -= 1
			j -= 1
			totAligned = totAligned + 1

        #Gap in S1, moving up
		if traceMat[i][j] == 1:
			s1_aligned = "_" + s1_aligned
			s2_aligned = s2[i-1][1] + s2_aligned
			i -= 1

        #Gap in S2, moving left
		if traceMat[i][j] == 2:
			s1_aligned = s1[j-1][1] + s1_aligned
			s2_aligned = "_" + s2_aligned
			j -= 1

	while i != 0 and j == 0:
		s1_aligned = "_" + s1_aligned
		s2_aligned = s2[i-1][1] + s2_aligned
		i -= 1

	while i == 0 and j != 0:
		s1_aligned = s1[j-1][1] + s1_aligned
		s2_aligned = "_" + s2_aligned
		j -= 1
    	
	return [s1_aligned, s2_aligned, totAligned]

#Temporal Smith Waterman score matrix calculation
#Here we fill our H matrix with values, starting from the top left, according to those equations defined in 10.1109/dsaa.2015.7344785
#We also store secondary sequences which achieved a high score, depending on the user's mem input, which can be used both for comparison
#and in the case that our shorter s1 sequence aligns with a high score at multiple, independent locations in s2
def TSW_scoreMat(s1,s1_len,s2,s2_len,g,T,H,TR,TC,traceMat,s,mem = 0):
	#initialise time penalty at 0
	tp = 0

	#initialise max_score and max_index
	max_score = -1
	max_index = (-1,-1)

	#Initialise memory
	mem_index = []
	mem_score = []

	for i in range(1,s2_len+1):
		for j in range(1,s1_len+1):
			if i == 1 and j == 1:
				tp = 0
			#Calculate correct time penalty according to ref.
			else:
				tpx = float(s2[i-1][0]) + float(TR[i-1][j-1])
				tpy = float(s1[j-1][0]) + float(TC[i-1][j-1])

				tpmax = max(tpx,tpy)

				#Avoid division by 0
				if tpmax == 0:
					tpmax = 1e-24
				
				tp = T*(abs(tpx-tpy)/tpmax)

			#Calculate potential Hmat choices
			score_match = score(s,s1[j-1][1],s2[i-1][1])
			Hup = H[i-1][j-1] + score_match - tp
			Hmid = H[i-1][j] - g
			Hbot = H[i][j-1] - g

			#Select absolute maximum of options
			#0 is STOP, 1 is DIAG, 2 is UP, 3 is LEFT
			H[i][j] = max([0,Hup,Hmid,Hbot])

			##Update traceMat
			#Tij = 0 for STOP
			#Tij = 1 for Match
			#Tij = 2 for Up
			#Tij = 3 for Left
			traceMat[i][j] = [0,Hup,Hmid,Hbot].index(H[i][j])

			traceVal = traceMat[i][j]
			matVal = H[i][j]

			#Update TR and TC matrices according to choice
			if matVal == 0:
				TR[i][j] = 0
				TC[i][j] = 0			
			elif traceVal == 1:
				TR[i][j] = 0
				TC[i][j] = 0
			elif traceVal == 2:
				TR[i][j] = TR[i-1][j] + float(s2[i-1][0])
				TC[i][j] = TC[i-1][j]
			elif traceVal == 3:
				TR[i][j] = TC[i][j-1]
				TC[i][j] = TC[i][j-1] + float(s1[j-1][0])

			if H[i][j] >= max_score:
				max_index = (i,j)
				max_score = H[i][j]

			#Add results to memory vector
			mem_index.append((i,j))
			mem_score.append(H[i][j])

	#Remove extra memory results
	mem_array = np.asarray(list(zip(mem_score,mem_index)), dtype=object)
	mem_array = mem_array[mem_array[:, 0].argsort()]

	#If mem == -1, we calculate the maximum number of possible matches when aligning s1 into s2
	#and keep all sequences with the corresponding score or higher
	if mem == -1:
		mem = max(1,math.floor(s2_len/s1_len))
		mem_min = mem_array[-mem][0]
		mem_array = mem_array[ mem_min <= mem_array[:,0] ]#
		mem_score, mem_index = mem_array[:, 0], mem_array[:, 1]

	#If mem == 0, we discard all memory indices and scores and return only the first encountered, highest scoring alignment
	elif mem == 0:
		mem_index = []
		mem_score = []

	#If mem >= 1, we return mem number of sequences
	else:
		mem_min = mem_array[-mem][0]
		mem_array = mem_array[ mem_min <= mem_array[:,0] ]
		mem_score, mem_index = mem_array[:, 0], mem_array[:, 1]

	#Reverse mem index such that the highest score/index is in position 0
	mem_score = mem_score[::-1]
	mem_index = mem_index[::-1]

	finalScore = max_score
	finalIndex = max_index
	return finalScore, finalIndex, mem_index, mem_score

#Here we utilise our traceback matrix to fill in all matches and gaps in each sequence to return the aligned sequences
#which achieved the highest scoring cell value in H
def align_TSW(traceMat, s1, s2, s1_len, s2_len, max_index):
	s1_aligned = ""
	s2_aligned = ""   
	temp_s1_aligned = ""   
	temp_s2_aligned = ""  
	max_j, max_i = max_index
	i = s1_len
	j = s2_len

	totAligned = 0

	# Tracing and computing the pathway with the local alignment
	while traceMat[max_j][max_i] > 0:
		#DIAG - match
		if traceMat[max_j][max_i] == 1:
			temp_s1_aligned = s1[max_i - 1][1] + s1[max_i - 1][0]
			temp_s2_aligned = s2[max_j - 1][1] + s2[max_j - 1][0]
			max_i -= 1
			max_j -= 1
			i -= 1
			j -= 1
			totAligned += 1
		
		#UP - gap in s2
		elif traceMat[max_j][max_i] == 3:
			temp_s1_aligned = s1[max_i - 1][1] + s1[max_i - 1][0]
			temp_s2_aligned = '__'
			max_i -= 1 
			i -= 1
		
		#LEFT - gap in s1
		elif traceMat[max_j][max_i] == 2:
			temp_s1_aligned = '__'
			temp_s2_aligned = s2[max_j - 1][1] + s2[max_j - 1][0]
			max_j -= 1
			j -= 1

		#Recombine temp sequences
		s1_aligned = s1_aligned + temp_s1_aligned
		s2_aligned = s2_aligned + temp_s2_aligned

	#Ensure that the final match character is included
	if max_i != 0 and max_j != 0:
		temp_s1_aligned = s1[max_i - 1][1] + s1[max_i - 1][0]
		temp_s2_aligned = s2[max_j - 1][1] + s2[max_j - 1][0]
		max_i -= 1
		max_j -= 1
		i -= 1
		j -= 1
		totAligned += 1

		s1_aligned = s1_aligned + temp_s1_aligned
		s2_aligned = s2_aligned + temp_s2_aligned

	# Reversing the order of the sequences
	s1_aligned = s1_aligned[::-1]
	s2_aligned = s2_aligned[::-1]

	return s1_aligned, s2_aligned, totAligned

def removeOverlaps():
	#Not implemented
	a = 0

#Here we perform all initialisation and alignment steps and output the corresponding data
def temporal_alignment(s1,s2,g,T,s,local_align,verbose,mem=-1,removeOverlap=0):
	#Initialise sequence lengths
	s1_len = len(s1)
	s2_len = len(s2)

	if s1_len > s2_len:
		print("Warning: Your regimen sequence appears to be longer than your patient's drug sequence. Consider checking patient record.")

	#Initialise matrices
	H = init_Hmat(s1_len,s2_len,g,local_align)
	TR = init_TRmat(s1,s1_len,s2,s2_len,local_align)
	TC = init_TCmat(s1,s1_len,s2,s2_len,local_align)
	traceMat = init_traceMat(s1_len,s2_len,local_align)
	secondary = 0

	#Perform global alignment
	if local_align == 0:
		if verbose == 1 or verbose == 2:
			print("Performing global alignment...")

		#Initialise return vars
		returnDat = [str(s1).strip('[]'),str(s2).strip('[]'),"NA","NA"]

		#Generate score matrix
		finalScore = TNW_scoreMat(s1,s1_len,s2,s2_len,g,T,H,TR,TC,traceMat,s)

		#Track alignments through the traceback matrix and append data to return DF
		s1_aligned, s2_aligned, totAligned = align_TNW(traceMat, s1, s2, s1_len, s2_len)
		returnDat.append([str(s1_aligned).strip('[]'),str(s2_aligned).strip('[]'),finalScore,totAligned])

		if verbose == 2:
			print("Final score matrix: ")
			print(H)
			print()
			print("Final traceback matrix: ")
			print(traceMat)
			print()

		print("Optimal alignment score: ")
		print(finalScore)
		print()
		print("Best global alignment of S1, S2")
		print(s1_aligned)
		print(s2_aligned)

		if removeOverlap == 1:
			print("Warning: You have attempted to remove overlapping alignments from a single global alignment, please check your settings.")

		return returnDat

	if local_align == 1:
		if verbose == 1 or verbose == 2:
			print("Performing local alignment...")

		#Initialise sequence length regex pattern
		pat = "[0-9][A-Z]|__"

		#Initialise return vars
		returnDat = [str(s1).strip('[]'),str(s2).strip('[]'),"","","","","","",""]
		returnDat = np.array(returnDat, dtype=object)

		#Generate score matrix
		finalScore, finalIndex, mem_index, mem_score = TSW_scoreMat(s1,s1_len,s2,s2_len,g,T,H,TR,TC,traceMat,s,mem)

		if verbose == 2:
			print("Final score matrix: ")
			H = pd.DataFrame(H)
			s1p = ["NA"] + [''.join(i) for i in s1]
			Hp = H.set_axis(s1p, axis = 1, copy=False)
			s2p = ["NA"] + [''.join(i) for i in s2]
			Hp = Hp.set_axis(s2p, axis = 0, copy=False)
			print(Hp)
			print()
			print("Final traceback matrix: ")
			traceMatp = pd.DataFrame(traceMat)
			traceMatp = traceMatp.set_axis(s1p, axis = 1, copy=False)
			traceMatp = traceMatp.set_axis(s2p, axis = 0, copy=False)
			print(traceMatp)
			print()

		#Check to see if any other alignments were stored in memory.
		#We need this extra if statement here to allow for the functionality of m = 0, which generates an empty mem_index
		if len(mem_score) > 1:
			#Track alignments through the traceback matrix
			s1_aligned, s2_aligned, totAligned = align_TSW(traceMat, s1, s2, s1_len, s2_len, mem_index[0])
			#Generate some summary data
			s_a_len = len(re.findall(pat,s1_aligned))
			s1_start = mem_index[0][1] - s_a_len
			s1_end = mem_index[0][1]
			s2_start = mem_index[0][0] - s_a_len
			s2_end = mem_index[0][0]

			#Append data to return DF
			returnDat = np.append(returnDat,([s1_aligned,s2_aligned,finalScore,s1_start,s1_end,
				s2_start,s2_end,s_a_len,totAligned]), axis = 0)
		else: 
			#Track alignments through the traceback matrix
			s1_aligned, s2_aligned, totAligned = align_TSW(traceMat, s1, s2, s1_len, s2_len, finalIndex)
			#Generate some summary data
			s_f_len = len(re.findall(pat,s1_aligned))
			s1_start = finalIndex[1] - s_f_len
			s1_end = finalIndex[1]
			s2_start = finalIndex[0] - s_f_len
			s2_end = finalIndex[0]

			#Append data to return DF
			returnDat = np.append(returnDat,([s1_aligned,s2_aligned,finalScore,s1_start,s1_end,
				s2_start,s2_end,s_f_len,totAligned]), axis = 0)

		if verbose == 1 or verbose == 2:
			print("Best local alignment of S1, S2")
			print(finalIndex)
			print(s1_aligned)
			print(s2_aligned)
			print("Score: ")
			print(finalScore)
			print()

		if len(mem_score) > 1:
			if verbose == 1 or verbose == 2:
				print("Secondary alignments:")
		
			for i in range(1,len(mem_index)):
				secondary = 1
				#Track alignments through the traceback matrix
				s1_aligned_t, s2_aligned_t, totAligned_t = align_TSW(traceMat, s1, s2, s1_len, s2_len, mem_index[i])
				#Generate some summary data
				s_a_t_len = len(re.findall(pat,s1_aligned_t))
				s1_start = mem_index[i][1] - s_a_t_len
				s1_end = mem_index[i][1]
				s2_start = mem_index[i][0] - s_a_t_len
				s2_end = mem_index[i][0]
				#Append data to return DF
				returnDat = np.append(returnDat,[s1_aligned_t,s2_aligned_t,mem_score[i],s1_start,s1_end,
					s2_start,s2_end,s_a_t_len,totAligned_t], axis = 0)	
	
				if verbose == 1 or verbose == 2:
					print(mem_index[i])
					print(s1_aligned_t)
					print(s2_aligned_t)
					print("Score: ")
					print(mem_score[i])
					print()

		#If we did have secondary alignments, reshape array accordingly
		if secondary == 1:
			returnDat = returnDat.reshape(len(mem_index)+1,9)
			returnDat = pd.DataFrame(returnDat)
		#If not, ensure return frame has the shape 2,9
		else:
			returnDat = returnDat.reshape(2,9)
			returnDat = pd.DataFrame(returnDat)

		#Remove overlaps
		if removeOverlap == 1:
			if verbose == 1 or verbose == 2:
				print("Removing overlaps...")

			rows, cols = np.shape(returnDat)
			#Get a list of starts and stops, ignoring row 1 (Containing the input)
			interval_List = returnDat[[5,6]].values.tolist()[1:]

			#Initialise two sets, one to keep track of rows to include and one to keep track of bases already covered by an interval
			keep_rows = [0]
			covered_bases = set([])
			#Initialise a row tracker
			i = 1

			#For each interval, check if the union of sets between the interval and the covered bases exists, if it does, we ignore 
			#that row, do not add it to our keep vector and do not add the bases it covers to the covered bases vector.
			for start, end in interval_List:
				interval = set(range(int(start),int(end)+1))

				if len(interval & covered_bases) == 0:
					keep_rows.append(i)
					covered_bases.update(interval)

				i += 1

			returnDat = returnDat.loc[keep_rows]

		return returnDat


