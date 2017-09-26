# Recombination rate estimation and recombination map reconstruction
Bioinformatics scripts for estimating per-generation recombination rate (r) from parent-child genotyping data.

The main idea is to employ the likelihood recursive function for kinship [1],  where r as a variable, and use Maximum Likelihood Estimation (MLE) to estimate r. Technologically, r between two SNPs can be inferred in condition that there is linkage disequilibrium (LD) between the two SNPs. However, over 15K parent-child pairs are usually needed to infer r reliably (sample size comparable with the deCODE project and the AA map project). Thus, to enhance the power, we usually include 10 serial SNPs (with a region of ~ 10kbp) for inferring the r at a 10kbp region.

We propose a Monte Carlo approach to infer the r between two SNPs with a evaluation of reliance. Tur script for converting the recombination rates with SNPs as the coordinate (gmap) to SRR(standard recombination rate) with physical distances as coordinate (rmap).

Usage:
We use the plink format *.ped , *.map and *.frq file in the following steps of estimation.

Plink format:

    -MRlt.ped:

        fam1"\t"id1"\t0\t0\t0\t0\tA T\tC G\n"
        fam1"\t"id2"\t0\tid1\t0\t0\tA T\tC G\n"
        fam2"\t"id3"\t0\t0\t0\t0\tA T\tC G\n"
        fam2"\t"id4"\t0\tid3\t0\t0\tA T\tC G\n"

    -Unr.ped:
        fam1"\t"id1"\t0\t0\t0\t0\tA T\tC G\n"
        fam2"\t"id3"\t0\t0\t0\t0\tA T\tC G\n"
        
    -*.map:
        chr"\t"rs"\t"0"\t"loc
        
    -*.frq:
        rs"\t"A"\t"T"\t"0.1577
        

Step1 (python):

Load the genotype dataset in a plink forma to estimate r between pairs of SNP sites by using MLE. The gametes frequencies of two SNPs are calucated by using the CubeX [2] python program. A total of 500 times of disturbance are performed for Monte Carlo estimation of r and reliance.

    Input parameters include: 
  
      -UnrPlinkPath(unrelated_plink_file),
  
      -ExrPlinkPath(related_plink_file)
  
      -OutPath(output)
    
      -skip_sites
  
      -ncpus(with pp module)

Step2 (r): 

We use the 'tmvtnorm' r-package[3] to estimate the mean r values of the 500 disturbation expriments, which form a truncated normal distribution. 

    Input parameters include:
    
      -input_path
      
      -output_path
    
step 3-5 (python): 

In these steps, we convert gmap to rmap as that in the deCode project[4]. We use the linear interpolation to calculate the SRR. 


    Input parameters include:
    
      -step3
      
      -skip1file: adjacent sites gmap file
      
      -skipnfile: when n, the file's result is calculated by the two sites whose interval has another n-1 site(s).
      
      -Peakpath : Repeat estimation at one site at one site.
      
      -pntspath : All distribution value at one site.
      
      
      -step4
      
      -inputpath : step3.py creates peakpath
      
      -outputpath: accurate gmap file


      -step5
      
      -inputPath : accurate gmap file
      
      -change in code (bin start, bin end, bin width)
      
      It will print the rmap value(cM/Mb)





1. Brenner, C. H. (1997). Symbolic Kinship Program. Genetics, 145(2), 535–542.

2. Gaunt, T. R., Rodríguez, S., & Day, I. N. (2007). Cubic exact solutions for the estimation of pairwise haplotype frequencies: implications for linkage disequilibrium analyses and a web tool “CubeX.” BMC Bioinformatics, 8(1), 428. https://doi.org/10.1186/1471-2105-8-428

3. Wilhelm, S., & Manjunath, B. G. (2010). tmvtnorm : A Package for the Truncated Multivariate Normal Distribution Generation of random numbers computation of marginal densities. The R Journal, 2(1), 25–29.

4. Kong, A., Gudbjartsson, D. F., Sainz, J., Jonsdottir, G. M., Gudjonsson, S. A., Richardsson, B., … Stefansson, K. (2002). 2002:A high-resolution recombination map of the human genome. Nature Genetics, 31(3), 241–7. https://doi.org/10.1038/ng917
