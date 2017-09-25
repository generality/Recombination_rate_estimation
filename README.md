# Recombination_rate_estimation
(Bioinformatics) Script for estimating per-generation RR(recombination rate).

Our main idea of estimating RR is MLE(Maximum Likelihood Estimation) via transforming the kinship test likelihood recursive function[1]. In the procedure of estiamting, we propose repeat disturbance experiment to curve the estiamtion's distribution and infer the parameter which is an approximately accurate estimation. Finally, we provide our script for converrting RR to SRR(standard recombination rate) to compare with other results.

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
        

Step1.py:

In step1, we use plink format genotype file to create a set of recombination rate estimation between pairs of SNP sites with MLE. One pair of sites has 500 times repeat disturbance experiment. What's more, we use 'CubeX' python script in MLE to help us calculate the pairs' gametes frequency[2]. 

    Input parameters include: 
  
      -UnrPlinkPath(unrelated_plink_file),
  
      -ExrPlinkPath(related_plink_file)
  
      -OutPath(output)
    
      -skip_sites
  
      -ncpus(with pp module)

Step2.r: 

We need a r-package named 'tmvtnorm'[3] to create gmap. In this step, we estimate the mean of disturbance experiment's data set whose distribution is truncated normal distribution. 

    Input parameters include:
    
      -input_path
      
      -output_path
    
Step3.py/Step4.py/Step5.py: 

In this step, we convert gmap to rmap. Like Decode genetic method[4], we use the linear interpolation to calculate the SRR. 


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
