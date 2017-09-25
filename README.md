# Recombination_rate_estimation
Bioinformatics_code

Before the step of estimation, you need to provide the plink format with tab and space.


Step1.py: We use the CubeX's code to calculate gametes' frequency[1]. 

    Input parameters include: 
  
      -UnrPlinkPath(unrelated_plink_file),
  
      -ExrPlinkPath(related_plink_file)
  
      -OutPath(output)
    
      -skip_sites
  
      -ncpus(with pp module)

Step2.r: We need a r-package named 'tmvtnorm'[2] to create gmap.

    Input parameters include:
    
      -input_path
      
      -output_path
    
Step3.py/Step4.py/Step5.py: convert gmap to rmap

    Input parameters include:
    
      -step3
      
      -skip1file: adjacent sites gmap file
      
      -skipnfile: when n, the file's result is calculated by the two sites whose interval has another n-1 site(s).
      
      -Peakpath : Repeat estimation at one site at one site.
      
      -pntspath : All distribution value at one site.
      
      
      -step4
      
      -inputpath : step3.py creates peakpath
      
      -outputpath: accurate gmap file


      -step5[3]
      
      -inputPath : accurate gmap file
      
      -change in code (bin start, bin end, bin width)
      
      It will print the rmap value(cM/Mb)








1. Gaunt, T. R., Rodríguez, S., & Day, I. N. (2007). Cubic exact solutions for the estimation of pairwise haplotype frequencies: implications for linkage disequilibrium analyses and a web tool “CubeX.” BMC Bioinformatics, 8(1), 428. https://doi.org/10.1186/1471-2105-8-428

2. Wilhelm, S., & Manjunath, B. G. (2010). tmvtnorm : A Package for the Truncated Multivariate Normal Distribution Generation of random numbers computation of marginal densities. The R Journal, 2(1), 25–29.

3. Kong, A., Gudbjartsson, D. F., Sainz, J., Jonsdottir, G. M., Gudjonsson, S. A., Richardsson, B., … Stefansson, K. (2002). 2002:A high-resolution recombination map of the human genome. Nature Genetics, 31(3), 241–7. https://doi.org/10.1038/ng917
