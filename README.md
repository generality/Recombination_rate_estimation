# Recombination_rate_estimation
Bioinformatics_code

Before the step of estimation, you need to provide the plink format with tab and space.


Step1.py: We use the CubeX's code to calculate gametes' frequency[]. 

  Input parameters includes: 
    -UnrPlinkPath(unrelated_plink_file),
    -ExrPlinkPath(related_plink_file)
    -OutPath(output)
    -skip_sites
    -ncpus(with pp module)
