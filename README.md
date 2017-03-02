PhD Data
============

My PhD data sets, code and experimentation results

The PhD is entitled *Interactive, incremental and multi-level exploration of large collections of images*.
It has been defended in November the 22nd, 2016 in the University Fran�ois-Rabelais of Tours.


I - Datasets
--------------
* Dinet_images/			Bentham word images and class ground-truth
* Wang_images/			Wang data set images (1000 images)
* breiman.txt			Breiman data set (5000 x 41 matrix) - class in the last column
* corel68k.txt			Corel68k data set (68,040 x 58 matrix) - class in the last column
* iris.txt			Iris data set (150 x 5 matrix) - class in the last column
* MIRFLICKR-25000_CLD.txt		MirFlickr-25k CLD data set (25,000 x 193 matrix) - imagename in the first column	
* NASA_APOD_CLD.txt		Nasa APoD CLD data set (6826 x 193 matrix) image date in the first column
* NGA_Images_CLD.txt		National Gallery of Art Open Access Images CLD data set (43,721 x 193 matrix) image name in the first column
* NGA_Images_EHD.txt		National Gallery of Art Open Access Images EHD data set (43,721 x 81 matrix) image name in the first column
* Wang_CLD.txt			Wang CLD data set ( 1000 x 193 matrix) 
* wdbc.txt			Wisconsin Diagnostic Breast Cancer data set (569 x 31 matrix)  - imagename in the last column


II - Graphs
--------------
1. Chapter_5_Experimentations/	
  | 
  |- _X_DATASET/		
		Contains graphs generated by (i) brute force algorithm, (ii) [Hacid et al, 2007] algorithm and (iii) the proposed approximate incremental algorithm. 
        For the proposed algorithms, results are given for neighbourhood order L={2,3,4}

2. Chapter_6_Experimentations/	
		Contains the the HRNG generated by our algorithm for five data sets. 

3. Chapter_11_Experimentations/	
		Contains the 12 graphs generated for the Chapter 11 use case 
		One graph for each (feature, distance) pair

4. star.json	Basic star graph with seven nodes and six edges


III - Code
--------------

1. RNG_python		[Python] Build a Relative Neighbourhood Graph (RNG) 
2. iRNG_Hacid		[C++] Incremantally build the RNG with Hacid et al. algorithm
3. iRNG_Exact		[C++] Incremantally build the RNG with the proposed exact algorithm
4. iRNG_Approximate	[C++] Incremantally build the RNG with the proposed approximate algorithm

5. BIRCH++			[C++] Incremantally build the proposed hiearchical and graph-based hybrid structure (HRNG)
6. LayoutHRNG		[C++] Layout each level of the HRNG with the Stress Majorization graph drawing algorithm
7. SplitHRNG		[C++] Split a huge inner HRNG file into several files, one for each internal nodes

8. sigma.layouts.maxentStress		[JavaScript] Plugins for the Linkurious JavaScript library
9. sigma.layouts.stressMajorization	[JavaScript] Plugins for the Linkurious JavaScript library
10. image_treemap					[JavaScript] Custom JavaScript object: Visual treemap
11. image_breadcrumb				[JavaScript] Custom JavaScript object: Visual breadcrumb
