import numpy as np
import pandas as pd
import math as math
import re

def lossFunction(T,t1,t2,method,i,j):
	absDiff = abs(t1-t2)
	maxT = max(t1,t2)

	#Avoid division by zero
	if maxT == 0:
		maxT = 1e-24

	#Ensure that trailing time penalty is set to an unrealistic minimum
	if j == 1:
		tp = 1e-24

	elif method == "PropDiff":
		if t1 == 0 or t2 == 0:
			tp = T * ((absDiff)/(maxT+1))

		else:
			tp = T * ((absDiff)/(maxT))

	elif method == "AbsDiff":
		tp = T * absDiff

	elif method == "Quadratic":
		tp = T * pow(absDiff,2)

	elif method == "PropQuadratic":
		tp = (T * pow(absDiff,2))/maxT

	elif method == "LogCosh":
		tp = math.cosh(math.log(absDiff))

	else:
		print("Bad method. Please choose one of `PropDiff`, `AbsDiff`, `Quadratic`, `PropQuadratic` or `LogCosh`")

	return tp


def score(s,x,y):
	score = s[x][y]

	return score

def swapPos(s2,i,k):

	temp = s2[i-1][1]

	s2[i-1][1] = s2[k-1][1]
	s2[k-1][1] = temp

	return s2

def TSW_scoreMat(s1,s1_len,s2,s2_len,g,T,H,TR,TC,traceMat,s,method):
	#initialise time penalty at 0
	tp = 0

	i = 1
	#Traverse drug record (drugRec=s2=i)
	while i < s2_len+1:
		j = 1
		#Traverse regimen (regimen=s1=j)
		while j < s1_len+1:
		
			#Dynamic re-ordering

			#Check that i does not refer to the start or end of the drug record
			if (i > 1 and i < s2_len):

				#Check that either record contains a 0 in the first slot
				if float(min(s1[j-1][0],s2[i-1][0])) == 0:

					#Check for a match between drug record and regimen
					if s1[j-1][1] == s2[i-1][1]:

						#Start of regimen case - no need to check for leading matches
						if j == 1 and float(s2[i-1][0]) == 0:

							k = i - 1
							s2 = swapPos(s2,i,k)

							i = i-1
							continue

						#End of regimen case - no need to check for trailing matches
						elif j == s1_len and float(s2[i-1][0]) == 0:

							k =  i + 1

							#Ensure that next entry into drug record is not a new day
							if float(s2[k-1][0]) == 0:

								s2 = swapPos(s2,i,k)

								i = i-1
								continue

			#Calculate time penalty
			if i == 1 and j == 1:
				tp = 0
			else:
				tpx = float(s2[i-1][0]) + float(TR[i-1][j-1])
				tpy = float(s1[j-1][0]) + float(TC[i-1][j-1])

				tp = lossFunction(T,tpx,tpy,method,i,j)				

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



			if j == s1_len and i < s2_len:
				if(float(s2[i][0]) < float(s1[0][0])):
					H[i][j] = max(H[i][j]-lossFunction(T,float(s2[i][0]),float(s1[0][0]),method,i,j),0)

			j += 1

		i += 1


def find_best_score(H,s2_len,s1_len,verbose):

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

	mem = max(1,s1_len - s2_len + 1)
	mem_min = mem_array[-mem][0] * 0.9
	mem_array = mem_array[ mem_min <= mem_array[:,0] ]
	mem_score, mem_index = mem_array[:, 0], mem_array[:, 1]

	if(verbose == 2):
		print("Calculated mem: ")
		print(mem)

	mem_score = mem_score[::-1]
	mem_index = mem_index[::-1]

	finalScore = max_score
	finalIndex = max_index

	return finalScore, finalIndex, mem_index, mem_score
