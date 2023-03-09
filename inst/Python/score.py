import numpy as np
import pandas as pd
import math as math
import re

def score(s,x,y):
	score = s[x][y]

	return score

def TSW_scoreMat(s1,s1_len,s2,s2_len,g,T,H,TR,TC,traceMat,s):
	#initialise time penalty at 0
	tp = 0

	for i in range(1,s2_len+1):
		for j in range(1,s1_len+1):

			#Calculate time penalty
			if i == 1 and j == 1:
				tp = 0
			else:
				tpx = float(s2[i-1][0]) + float(TR[i-1][j-1])
				tpy = float(s1[j-1][0]) + float(TC[i-1][j-1])

				tpmax = max(tpx,tpy)

				#Avoid division by zero
				if tpmax == 0:
					tpmax = 1e-24
				
				tp = T*(abs(tpx-tpy)/tpmax)

				if j == 1:
					tp = 1e-24

			#Calculate score
			score_match = score(s,s1[j-1][1],s2[i-1][1])

			#Calculate potential values of H ij via score, tp and Hi-1, Hj-1
			Hup = H[i-1][j-1] + score_match - tp
			Hmid = H[i-1][j] - g
			Hbot = H[i][j-1] - g

			#Select move
			H[i][j] = max([0,Hup,Hmid,Hbot])

			#Record move
			traceMat[i][j] = [0,Hup,Hmid,Hbot].index(H[i][j])

			traceVal = traceMat[i][j]
			matVal = H[i][j]

			#ZERO
			#if matVal == 0:
			#	TR[i][j] = 0
			#	TC[i][j] = 0	

			#DIAGONAL		
			if traceVal == 1:
				TR[i][j] = 0
				TC[i][j] = 0

			#HORIZONTAL
			elif traceVal == 2:
				TR[i][j] = TR[i-1][j] + float(s2[i-1][0])
				TC[i][j] = TC[i-1][j]

			#VERTICAL
			elif traceVal == 3:
				TR[i][j] = TC[i][j-1]
				TC[i][j] = TC[i][j-1] + float(s1[j-1][0])


def find_best_score(H,s2_len,s1_len,mem,verbose):
	mem_index = []
	mem_score = []

	max_score = -1
	max_index = (-1,-1)

	for i in range(0,s1_len+1):
		mem_index.append((i,s2_len))
		mem_score.append(H[i][s2_len])

		if H[i][s2_len] >= max_score:
				max_index = (i,s2_len)
				max_score = H[i][s2_len]

	mem_array = np.asarray(list(zip(mem_score,mem_index)), dtype=object)
	mem_array = mem_array[mem_array[:, 0].argsort()]

	#Select from mem according to the requested number of alignments
	#Mem = -1 : As many non-overlapping alignments as possible
	if mem == -1:
		mem = max(1,math.ceil(s1_len/s2_len))
		mem_min = mem_array[-mem][0]
		mem_array = mem_array[ mem_min <= mem_array[:,0] ]#
		mem_score, mem_index = mem_array[:, 0], mem_array[:, 1]
		if(verbose == 2):
			print("Calculated mem: ")
			print(mem)

	#Mem = 0 : Exactly 1 alignment
	elif mem == 0:
		mem_index = []
		mem_score = []

	# Mem >= 1 : Return N alignments and all alignments with the same score as the Nth alignment 
	else:
		mem_min = mem_array[-mem][0]
		mem_array = mem_array[ mem_min <= mem_array[:,0] ]
		mem_score, mem_index = mem_array[:, 0], mem_array[:, 1]

	mem_score = mem_score[::-1]
	mem_index = mem_index[::-1]

	finalScore = max_score
	finalIndex = max_index

	return finalScore, finalIndex, mem_index, mem_score
