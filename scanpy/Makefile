

15days.h5ad: 
	python preprocess.py samples15d.tsv 15days
	python split.py 15days.h5ad 

35days.h5ad: 
	python preprocess.py samples35d.tsv 35days 

photoreceptors.h5ad: 
	python merge.py 15days.h5ad 35days.h5ad photoreceptors.h5ad
