import numpy as np
import pandas as pd
import math as math
import re

def score(s,x,y):
	score = s[x][y]

	return score

def TNW_scoreMat(s1,s1_len,s2,s2_len,g,T,H,TR,TC,traceMat,s):
	tp = 0
	finalScore = 0

	for i in range(1,s2_len+1):
		for j in range(1,s1_len+1):

			if i == 1 and j == 1:
				tp = 0
			elif j == 1:
				print("J reset")
				tp = 0
			elif i == 1:
				print("I reset")
			else:
				tpx = float(s2[i-1][0]) + float(TR[i-1][j-1])
				tpy = float(s1[j-1][0]) + float(TC[i-1][j-1])

				tpmax = max(tpx,tpy)

				if tpmax == 0:
					tpmax = 1e-8
				
				tp = T*(abs(tpx-tpy)/tpmax)

			score_match = score(s,s1[j-1][1],s2[i-1][1])
			Hup = H[i-1][j-1] + score_match - tp
			Hmid = H[i-1][j] - g
			Hbot = H[i][j-1] - g

			H[i][j] = max([Hup,Hmid,Hbot])
			finalScore += H[i][j]

			traceMat[i][j] = [Hup,Hmid,Hbot].index(H[i][j])
			traceVal = traceMat[i][j]

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

def TSW_scoreMat(s1,s1_len,s2,s2_len,g,T,H,TR,TC,traceMat,s,mem = 0):
	#initialise time penalty at 0
	tp = 0
	max_score = -1
	max_index = (-1,-1)

	mem_index = []
	mem_score = []

	for i in range(1,s2_len+1):
		for j in range(1,s1_len+1):

			if i == 1 and j == 1:
				tp = 0
			else:
				tpx = float(s2[i-1][0]) + float(TR[i-1][j-1])
				tpy = float(s1[j-1][0]) + float(TC[i-1][j-1])

				tpmax = max(tpx,tpy)

				if tpmax == 0:
					tpmax = 1e-24
				
				tp = T*(abs(tpx-tpy)/tpmax)

				if j == 1:
					tp = 1e-24

			score_match = score(s,s1[j-1][1],s2[i-1][1])
			Hup = H[i-1][j-1] + score_match - tp
			Hmid = H[i-1][j] - g
			Hbot = H[i][j-1] - g

			H[i][j] = max([0,Hup,Hmid,Hbot])

			traceMat[i][j] = [0,Hup,Hmid,Hbot].index(H[i][j])

			traceVal = traceMat[i][j]
			matVal = H[i][j]

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

			mem_index.append((i,j))
			mem_score.append(H[i][j])

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