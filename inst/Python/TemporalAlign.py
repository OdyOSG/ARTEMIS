import numpy as np
import pandas as pd
import math as math

#Initialise H matrix that will be used for scoring according to sequence size
def init_Hmat(s1_len,s2_len,g,local_align):

	if local_align == 0:
		H = np.zeros((s2_len+1,s1_len+1), float)

		#Fill the first column with incrementing gap penalties
		for row in range(s2_len+1):
			H[row][0] = -row*g
	
		#Fill the first row with incrementing gap penalties
		for col in range(s1_len+1):
			H[0][col] = -col*g
	
		return H

	elif local_align == 1:
		H = np.zeros((s2_len+1,s1_len+1), float)

		#Fill the first column with incrementing gap penalties
		for row in range(s2_len+1):
			H[row][0] = 0
	
		#Fill the first row with incrementing gap penalties
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

def init_traceMat(s1_len,s2_len,local_align):
	if local_align == 0:
		traceMat = np.zeros((s2_len+1,s1_len+1), float)
		#Fill top row with 2
		traceMat[0] = 2
		#Fill left col with 1
		traceMat[:,0] = 1
		#Fill in DONE square
		traceMat[0][0] = -1
		
	elif local_align == 1:
		traceMat = np.zeros((s2_len+1,s1_len+1), float)
		#Fill top row with 0
		traceMat[0] = 0
		#Fill left col with 0
		traceMat[:,0] = 0
		#Fill in DONE square
		traceMat[0][0] = 0

	return traceMat

def score(s,x,y):
	score = s[x][y]

	return score

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

	if mem == -1:
		mem = max(1,math.floor(s2_len/s1_len))
		mem_min = mem_array[-mem][0]
		mem_array = mem_array[ mem_min <= mem_array[:,0] ]#
		mem_score, mem_index = mem_array[:, 0], mem_array[:, 1]

	elif mem == 0:
		mem_index = []
		mem_score = []

	else:
		mem_min = mem_array[-mem][0]
		mem_array = mem_array[ mem_min <= mem_array[:,0] ]
		mem_score, mem_index = mem_array[:, 0], mem_array[:, 1]

	mem_score = mem_score[::-1]
	mem_index = mem_index[::-1]

	finalScore = max_score
	finalIndex = max_index
	return finalScore, finalIndex, mem_index, mem_score

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

def temporal_alignment(s1,s2,g,T,s,local_align,verbose,mem=-1):
	#Initialise sequence lengths
	s1_len = len(s1)
	s2_len = len(s2)

	#Initialise matrices
	H = init_Hmat(s1_len,s2_len,g,local_align)
	TR = init_TRmat(s1,s1_len,s2,s2_len,local_align)
	TC = init_TCmat(s1,s1_len,s2,s2_len,local_align)
	traceMat = init_traceMat(s1_len,s2_len,local_align)
	secondary = 0

	#Perform global alignment
	if local_align == 0:

		#Initialise return vars
		returnDat = [str(s1).strip('[]'),str(s2).strip('[]'),"NA","NA"]

		finalScore = TNW_scoreMat(s1,s1_len,s2,s2_len,g,T,H,TR,TC,traceMat,s)
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

		return returnDat

	if local_align == 1:

		#Initialise return vars
		returnDat = [str(s1).strip('[]'),str(s2).strip('[]'),"NA","NA","NA"]
		returnDat = np.array(returnDat, dtype=object)

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
			print(traceMat)
			print()

		if len(mem_score) > 1:
			s1_aligned, s2_aligned, totAligned = align_TSW(traceMat, s1, s2, s1_len, s2_len, mem_index[0])
			returnDat = np.append(returnDat,([s1_aligned,s2_aligned,finalScore,str(mem_index[0]),totAligned]), axis = 0)
		else: 
			s1_aligned, s2_aligned, totAligned = align_TSW(traceMat, s1, s2, s1_len, s2_len, finalIndex)
			returnDat = np.append(returnDat,([s1_aligned,s2_aligned,finalScore,str(finalIndex),totAligned]), axis = 0)

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
				s1_aligned_t, s2_aligned_t, totAligned_t = align_TSW(traceMat, s1, s2, s1_len, s2_len, mem_index[i])
	
				returnDat = np.append(returnDat,[s1_aligned_t,s2_aligned_t,mem_score[i],str(mem_index[i]),totAligned_t], axis = 0)	
	
				if verbose == 1 or verbose == 2:
					print(mem_index[i])
					print(s1_aligned_t)
					print(s2_aligned_t)
					print("Score: ")
					print(mem_score[i])
					print()



		if secondary == 1:
			returnDat = returnDat.reshape(len(mem_index)+1,5)
		else:
			returnDat = returnDat.reshape(2,5)


		


		return returnDat


