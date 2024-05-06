from re import sub

def align_TSW(traceMat, s1, s2, s1_len, s2_len, max_index):
	s1_aligned = ""
	s2_aligned = ""   
	temp_s1_aligned = ""   
	temp_s2_aligned = ""  
	max_j, max_i = max_index
	i = s1_len
	j = s2_len

	totAligned = 0

	while traceMat[max_j][max_i] > 0:
		if traceMat[max_j][max_i] == 1:
			temp_s1_aligned = s1[max_i - 1][1] + "." + s1[max_i - 1][0] + ";"
			temp_s2_aligned = s2[max_j - 1][1] + "." + s2[max_j - 1][0] + ";"
			max_i -= 1
			max_j -= 1
			i -= 1
			j -= 1
			totAligned += 1
		
		elif traceMat[max_j][max_i] == 3:
			temp_s1_aligned = s1[max_i - 1][1] + "." + s1[max_i - 1][0] + ";"
			temp_s2_aligned = '__' + ";"
			max_i -= 1 
			i -= 1
		
		elif traceMat[max_j][max_i] == 2:
			temp_s1_aligned = '__' + ";"
			temp_s2_aligned = s2[max_j - 1][1] + "." + s2[max_j - 1][0] + ";"
			max_j -= 1
			j -= 1

		s1_aligned = s1_aligned + temp_s1_aligned
		s2_aligned = s2_aligned + temp_s2_aligned

	if max_i != 0 and max_j != 0:
		temp_s1_aligned = s1[max_i - 1][1] + "." + s1[max_i - 1][0] + ";"
		temp_s2_aligned = s2[max_j - 1][1] + "." + s2[max_j - 1][0] + ";"
		max_i -= 1
		max_j -= 1
		i -= 1
		j -= 1
		totAligned += 1

		s1_aligned = s1_aligned + temp_s1_aligned
		s2_aligned = s2_aligned + temp_s2_aligned

	s1_list = s1_aligned.split(";")
	s2_list = s2_aligned.split(";")	

	s1_list_rev = [sub(r"([aA-zZ|~]+)\.([0-9]+)",r"\2.\1",i) for i in s1_list]
	s2_list_rev = [sub(r"([aA-zZ|~]+)\.([0-9]+)",r"\2.\1",i) for i in s2_list]

	s1_aligned = ";".join(s1_list_rev[::-1])
	s2_aligned = ";".join(s2_list_rev[::-1])

	return s1_aligned, s2_aligned, totAligned