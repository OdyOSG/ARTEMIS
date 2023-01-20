import numpy as np
import pandas as pd
import math as math
import re

pd.options.display.max_columns = None

def temporal_alignment(s1,s2,g,T,s,local_align,verbose,mem=-1,removeOverlap=0):
	s1_len = len(s1)
	s2_len = len(s2)
	if s1_len > s2_len:
		print("Warning: Your regimen sequence appears to be longer than your patient's drug sequence. Consider checking patient record.")

	H = init_Hmat(s1_len,s2_len,g,local_align)
	TR = init_TRmat(s1,s1_len,s2,s2_len,local_align)
	TC = init_TCmat(s1,s1_len,s2,s2_len,local_align)
	traceMat = init_traceMat(s1_len,s2_len,local_align)
	secondary = 0

	if local_align == 0:
		if verbose == 1 or verbose == 2:
			print("Performing global alignment...")

		returnDat = [str(s1).strip('[]'),str(s2).strip('[]'),"NA","NA"]
		finalScore = TNW_scoreMat(s1,s1_len,s2,s2_len,g,T,H,TR,TC,traceMat,s)
		s1_aligned, s2_aligned, totAligned = align_TNW(traceMat, s1, s2, s1_len, s2_len)
		returnDat.append([str(s1_aligned).strip('[]'),str(s2_aligned).strip('[]'),finalScore,totAligned])
		
		if verbose == 2:
			print("Final score matrix: ")
			print(H)
			print("Final traceback matrix: ")
			print(traceMat)

		print("Optimal alignment score: ")
		print(finalScore)
		print("Best global alignment of S1, S2")
		print(s1_aligned)
		print(s2_aligned)
		
		if removeOverlap == 1:
			print("Warning: You have attempted to remove overlapping alignments from a single global alignment, please check your settings.")

		return returnDat

	if local_align == 1:
		if verbose == 1 or verbose == 2:
			print("Performing local alignment...")

		pat = "[0-9][A-Z]|__"
		returnDat = [str(s1).strip('[]'),str(s2).strip('[]'),"","","","","","",""]
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
			print("Final traceback matrix: ")
			traceMatp = pd.DataFrame(traceMat)
			traceMatp = traceMatp.set_axis(s1p, axis = 1, copy=False)
			traceMatp = traceMatp.set_axis(s2p, axis = 0, copy=False)
			print(traceMatp)

		if len(mem_score) > 1:
			s1_aligned, s2_aligned, totAligned = align_TSW(traceMat, s1, s2, s1_len, s2_len, mem_index[0])
			s_a_len = len(re.findall(pat,s1_aligned))
			s1_start = mem_index[0][1] - s_a_len
			s1_end = mem_index[0][1]
			s2_start = mem_index[0][0] - s_a_len
			s2_end = mem_index[0][0]

			returnDat = np.append(returnDat,([s1_aligned,s2_aligned,finalScore,s1_start,s1_end,
				s2_start,s2_end,s_a_len,totAligned]), axis = 0)
		else: 
			s1_aligned, s2_aligned, totAligned = align_TSW(traceMat, s1, s2, s1_len, s2_len, finalIndex)
			s_f_len = len(re.findall(pat,s1_aligned))
			s1_start = finalIndex[1] - s_f_len
			s1_end = finalIndex[1]
			s2_start = finalIndex[0] - s_f_len
			s2_end = finalIndex[0]

			returnDat = np.append(returnDat,([s1_aligned,s2_aligned,finalScore,s1_start,s1_end,
				s2_start,s2_end,s_f_len,totAligned]), axis = 0)

		if verbose == 1 or verbose == 2:
			print("Best local alignment of S1, S2")
			print(finalIndex)
			print(s1_aligned)
			print(s2_aligned)
			print("Score: ")
			print(finalScore)

		if len(mem_score) > 1:
			if verbose == 1 or verbose == 2:
				print("Secondary alignments:")
		
			for i in range(1,len(mem_index)):
				secondary = 1
				s1_aligned_t, s2_aligned_t, totAligned_t = align_TSW(traceMat, s1, s2, s1_len, s2_len, mem_index[i])
				s_a_t_len = len(re.findall(pat,s1_aligned_t))
				s1_start = mem_index[i][1] - s_a_t_len
				s1_end = mem_index[i][1]
				s2_start = mem_index[i][0] - s_a_t_len
				s2_end = mem_index[i][0]
				returnDat = np.append(returnDat,[s1_aligned_t,s2_aligned_t,mem_score[i],s1_start,s1_end,
					s2_start,s2_end,s_a_t_len,totAligned_t], axis = 0)	
	
				if verbose == 1 or verbose == 2:
					print(mem_index[i])
					print(s1_aligned_t)
					print(s2_aligned_t)
					print("Score: ")
					print(mem_score[i])
					print()

		if secondary == 1:
			returnDat = returnDat.reshape(len(mem_index)+1,9)
			returnDat = pd.DataFrame(returnDat)
		else:
			returnDat = returnDat.reshape(2,9)
			returnDat = pd.DataFrame(returnDat)

		if removeOverlap == 1:
			if verbose == 1 or verbose == 2:
				print("Removing overlaps...")

			rows, cols = np.shape(returnDat)
			interval_List = returnDat[[5,6]].values.tolist()[1:]

			keep_rows = [0]
			covered_bases = set([])
			i = 1

			for start, end in interval_List:
				interval = set(range(int(start),int(end)+1))

				if len(interval & covered_bases) == 0:
					keep_rows.append(i)
					covered_bases.update(interval)

				i += 1

			returnDat = returnDat.loc[keep_rows]

		if verbose == 1 or verbose == 2:
			print("Done")
		
		return returnDat