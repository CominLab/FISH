#FISH - Fast Indexing for Spaced-seed Hashing

##Abstract

###Background  
Spaced-seeds, i.e. patterns in which some fixed positions are allowed to be wild-cards, play a crucial role in several bioinformatics applications involving substrings counting and indexing, by often providing better sensitivity with respect to k-mers based approaches. K-mers based approaches are usually fast, being based on efficient hashing and indexing that exploits the large overlap between consecutive k-mers. Spaced-seeds hashing is not as straightforward and it is usually computed from scratch for each position in the input sequence.
Recently, the FSH (Fast Spaced seed Hashing) approach was proposed to improve the time required for computation of the spaced seed hashing of DNA sequences with a speed-up of about 1.5 with respect to standard hashing computation.  

###Results
In this work we propose a novel algorithm, Fast Indexing for Spaced seed Hashing (FISH), based on the indexing of small blocks that can be combined to obtain the hashing of spaced-seeds of any length. The method exploits the fast computation of the hashing of runs of consecutive 1 in the spaced seeds, that basically correspond to k-mer of the length of the run. 

###Conclusions
We run several experiments on NGS data from metagenomic experiments, to assess the time required for the computation of the hashing for each position in each read with respect to several spaced seeds. In our experiments, FISH can compute the hashing values of spaced seeds with a speedup, with respect to the traditional approach, between 1.9x to 6.4x, depending on the structure of the spaced seed.  


#Download
*NB:This software tests only the computation time and give the results in files*  
[FISH_download](https://bitbucket.org/samu661/fish/downloads/FISH.tar.gz)  

##Compilation:
Open terminal and go to FISH/Release/ and then:  
make all  
  
If you need to test in parallel mode you must add -fopenmp option on compilation in every makefiles in Release's subdirectories.

###Option compilation
Y = Number of core  
-jY

#Algorithm Option
In the following paragraph is described the input file's structure and the parameters available in the FISH algorithm.  

##File accepted and structure
File accepted have the following structure:  
Structure file .fna example:  
> \>IDENTIFICATION  
> ATAATTGGCAAGTGTTTTAGTCTTAGAGAGATTCTCTAAGTCTAACTTGAACTCAATTTGGAATCATTTCCCAATTTTTA

Structure .fastq example:  
> @IDENTIFICATION #0/1  
> CCCATGCCTTTAGCCAAATTCACGGTTTGATCACCCCTAAAACCAGCCAATATACCGAAGTGGAAGCCAGCATAAATGGCCTCAATATTACCGAAATGGAT  
> +  
> HBIIIIIIIHHDIHIGIIGGIHIIGIDIIIIBIHI@IIH@HIIHIIF5IIHEII>BDAHIBIEDBEIDG@HAEH*I@AEI=#CE?G17EEDHDEB@@?#8B  
  
In this VERSION, paired-end reads are passed to the algorithm in separeted file in which the reads are paired 
in the same order in which are writen. So we raccomend to control the paired-end read if they are paired in the
correct manner.

##Parameter
**-si** File path single-end reads  
**-pi** File paths paired-end reads   
**-dirO** Path directory output files. Default: output/  
**-q** Enter a spaced seeds path as -q <AbsPathFile>. Every file's line must contain a spaced seed. Ex. 1\*\*\*1\*111. 1 is the simbol considered, any others are not valid.  
Default spaced seeds are:  
1111011101110010111001011011111 -> CLARK-S paper  
1111101011100101101110011011111 -> CLARK-S paper  
1111101001110101101100111011111 -> CLARK-S paper  
1111010111010011001110111110111 -> rasbhari minimizing overlap complexity  
1110111011101111010010110011111 -> rasbhari minimizing overlap complexity  
1111101001011100111110101101111 -> rasbhari minimizing overlap complexity  
1111011110011010111110101011011 -> rasbhari maximizing sensitivity  
1110101011101100110100111111111 -> rasbhari maximizing sensitivity  
1111110101101011100111011001111 -> rasbhari maximizing sensitivity  

#Run
Calls algorithm where is compiled:  
./FISH -si ../TestInputFile/long_example_1.fna -q 11101110110110111101  
./FISH -pi ../TestInputFile/short_example_1.fna.1 ../TestInputFile/short_example_1.fna.2 -q 11101110110110111101  
./FISH -pi ../TestInputFile/short_example_2.fna.1 ../TestInputFile/short_example_2.fna.2 -q 11101110110110111101 .dirO /home/user/desktop/  

#Publication
S.Girotto, M.Comin, C.Pizzi  
*Efficient computation of spaced seed hashing with block indexing*  
BMC Bioinformatics 2018, 19 (Suppl 15) :441  
DOI: https://doi.org/10.1186/s12859-018-2415-8  
  
#License
[MIT](https://bitbucket.org/samu661/fish/src/master/license.md)