# GB-Hyperlinks_N-S-Divide
Data and Code Repository to accompany the paper: 'Assessing the North-South Divide in Britain’s Digital Landscape through Stable Dynamic Embedding of Spatial Web Data'. 

The datasets used are:
*the Geoindex and Host Link Graphs which can be downloaded from JISC UK Web Domain Dataset: `https://data.webarchive.org.uk/opendata/ukwa.ds.2/`. This repository contains the code used to process the files, but not the files themselves. 
* ONS postcode lookup, which has many versions that can be downloaded from: https://geoportal.statistics.gov.uk/search?q=postcode%20lookup. (This project accessed a Nov2021 version)
* ONS file which tells you which Local Authority District (LAD) and Region each Middle Layer Super Output Area (MSOA) belongs to. This was accessed on 22-05-2022 at https://geoportal.statistics.gov.uk/datasets/national-statistics-postcode-lookup-november-2021/about. 


**01Bash:**
* ```MSOA_REG11NM_lookup_bash.sh``` outputs a file called ```MSOA_REG11NM_lookup.csv``` which has a column for MSOA and Region. It uses the file `Output_Area_to_LSOA_to_MSOA_to_Local_Authority_District_(December_2017)_Lookup_with_Area_Classifications_in_Great_Britain.csv` downloaded from: https://geoportal.statistics.gov.uk/datasets/national-statistics-postcode-lookup-november-2021/about
* The same file from the ONS is used to make a postcode (pcds) to MSOA (msoa11) lookup, and is done via the script `MSOA_all_nov2021_bash.sh` and the output is called `MSOA_UK_nov2021.csv`.
* `filter1PC.R` takes in the *.tsv* files downloaded from the Geoindex files at https://data.webarchive.org.uk/opendata/ukwa.ds.2/ and filters them for hosts with one postcode associated and saves as `all.1.csv`. 
* `lookup_bash.sh` takes in `all.1.csv` to get out the host and pc only called `lookup_all.csv` -- `all.1.csv` is a file from Emmanouil's previous project of all archived websites with one postcode. 

**02DataProcessing**: Script to make lookup table (host/postcode)
* `lookup_all_script.R` uses output `host-pc.csv` to make lookup table called `lookup_all.csv`. Keeps only one unique pc for a host.

**04HPCScripts:** 
* `MSOA_1996-2003_2005-2010_filter_couk_couk.sh` runs the R script called `1996-2003_2005-2010_couk_couk_filter.R` to filter for each year in the linkage data, links between 2 ".co.uk" hosts. THe output files take the form `<year>-couk-couk-linkage.tsv`.
* * ```MSOA_2005-2010_make_adjacencyMatrices_v2.sh``` uses the R script ```MSOA_2005-2010_MakeAdjMatr_v2.R``` , which has 2 different filters to remove self-hyperlinks (within the website). It does this by comparing host1 and host2. It outputs A_<year>_v2.mtx files for 2005-2010 inclusive. Changes made the filters in this document on 22/04/2022. This takes the `<year>-couk-couk-linkage.tsv` files and the MSOA-postcode lookup and matches both sides of the *'.co.uk'* to *'.co.uk'* links to an MSOA via their postcode, then aggregates this data by MSOA to output matrices `A_for_<year>_v2.mtx` files, where $A_{ij} represents the number of hyperlinks from MSOA[i] to MSOA[j] -- so these are weighted, directed matrices. Running script also prints an output of general summary of the A matrices to the console. These matrices represent **how many hyperlinks between MSOAs** - see **folder 12** for what we go on to use. 
  
**12BinaryHyperlinkConnections:**
* Here, instead of aggregating **how many hyperlinks between places**, we are making matrices of **how many website to website connections between places** are there. This means replacing the *# hyperlinks* by a 1 if non-zero and a 0 if zero, to make an adjacency matrix. Can think of this as taking the huge adj. matrix of website level and replacing all non-zero values with a 1. Then they are aggregated to MSOA and then LAD level by summing. `MSOA_2005-2010_makeAdjMatr_binary.R` does this. 
* Then there is a Python script `makeA_LAD_binary.ipynb`, that aggregates the MSOA level matrix to LAD level for just England and Wales LADs. (Also Winsorizing the matrices and log10 on them is done).
* The data to go from MSOA to LAD (for E&W) comes from: https://geoportal.statistics.gov.uk/datasets/middle-layer-super-output-area-2001-to-middle-layer-super-output-area-2011-to-local-authority-district-2011-lookup-in-england-and-wales/explore and is saved as `MSOA2001_MSOA2011_LAD2011_Lookup_EW.csv`. These files are often updated and renewer version uploaded, but searching "pcd msoa lsoa lad uk" should return the most similar data. 

**14NorthSouthDivide:**
* **LADS_list.csv** tells you the order of the indexes of the LAD x LAD matrices (the ones created via the code in **folder 12**). 
* `simulations.ipynb` runs a model with a SBM of 2 communities. It shows how the permutation tests look when you know what the outcome should look like. This helps to justify the plots, add explainability when we apply these methods to real data (where we do not know the underlying distribution that the nodes come from). The pdf plots in this folder are those generated in this notebook and are found in the Supplementary material Simulated Example.
* `NorthSouthDivide.ipynb` first contains the code to map the LAD x LAD matrices to Reg x Reg (Reg for region) matrices for E&W. These Region and LAD level matrices can be found in the folder too. There is some additional code that tests code from `simulations.ipynb` on the real data. 
  
**15 UKincScot:**
Note: I've called it UK everywhere, but it's actually GB as NI is not included. v5 is binary still, just includes GB. 

* `makeA_LAD_v5_UKincScot.ipynb` takes in the file `PCD_OA_LSOA_MSOA_LAD_FEB20_UK_LU.csv` which was downloaded from a link given in folder 12. The start of the notebook makes this file into a MSOA to LAD lookup table for the UK, which it saves as the file `MSOA_LAD_UK_lookup.csv` which has been saved and then can just be used (this is where it is made). The notebook makes the LADxLAD matrices (v5) for the UK, so Eng, Wales + Scot, using the v4 binary MSOA x MSOA matrices made in **folder 12**.
* `LADS_list_UK.csv` is a list of all the LADs in E, W & S. `A_LAD_<year>_UK_v5.mtx` are the matrices for the whole UK at LAD x LAD level. These are relabelled in the next folder and called `A_LAD_<year>_GB.mtx` to avoid confusion from mislabelling as UK when it is truly GB. `A_LAD_<year>_GB.mtx` is the name of these matrices which is used in all later folders. 
* `matrixSVD_plots_UK_binaryLog10Wins.R` uses the shp file from: 
https://geoportal.statistics.gov.uk/datasets/ons::local-authority-districts-december-2016-generalised-clipped-boundaries-in-the-uk-1/explore?location=53.588449%2C-0.820715%2C5.52
and does the SVD on the v5 A LADxLAD matrices. Makes the "animationData" and saves this as `animationDataUK.csv`. 
Second part of the `R` script makes plots: scree plots, 2D reference map plots, embedding plots. 
Note: the below file doesn't transform the A matrices, so is not the results reported in the paper, but is kept as it outputs a file used in **folder 18**
* `permutationTestsUK.ipynb` script to run permuation tests for GB. (Also has code to output a `LAD_NorthIndicator_GB_lookup.csv` file, which is copied into **folder 18** and used there)
  
**15-1 WinsorizationPoint:**
* `Winsorizing_beforeAfter_concatenation.R` implements a distance measure for the embeddings and compares for concatenating then Winsorizing and vice versa.
  
**16GB_ConcThenWins:**
* `matrixSVD_plots_UK_binaryLog10ConcWins.R` file takes in the LAD level A matrices, then concatenates, logs, then Winsorizes. Then does SVD to calculate embedding, and makes some plots of the embedding and a scree plot, and the file `animationDataGB_ConcWins.csv`.
* `permutationTestsGB_ConcWins.ipynb` uses the animation data to do the permutation tests for N/S, U/R, L/nL and plots graphs and outputs the values of the p-values. 
* `transformationsAPlots_GB_ConcWins.R` generates  plots to show the effects of concatenating then Winsorizing the data for GB.
  
**18Procrustes:**
* `ProcrustesGB_X_avg_v2.R` has lots of work in it, makes the plots and does the iterative Procrustes transformation. Does the zoom in on London, as well as doing some X_early and X_late Procrustes transformations. - Think I need to extract this information to then be able to apply the permutation tests on the data. 
* `Procust_early_late_for_PermTests.R` is a script to procrust the embedding into $C = (long, lat,  no.enterprises, zeros)$. Method (a) procrustes $X_{all average}$ into $C$ to learn the transformation parameters $s, R, tt$. It then uses these values to Procrust $X_{early}$ and $X_{late}$ into $C$. Method (b) does: $X_{early} \in \mathbb{R}^{n \times 4}$ and $X_{late}$ into $C = (long, lat,  no.enterprises, zeros)$. Then export the transformed dataframes as csv files to do some permutation tests on. The output csv files are copied into **20** to run permutation tests on them.
* `CheckNumberOfLinksPerYear_12Sep2023.ipynb` - plot to show that the number of links is the pattern we see in dim1 of the UASE embedding of the data (Figure 5 in the paper)
* `Procrust_NS_divide_plots.R` is done to make the plots for the paper wrt the Procrusted Dimensions rather than the embedding dimensions for better reading. `YhatProcrusted_allYears.csv` contains the Procrusted embeddings for each year all stacked into one matrix. 
* `Procrusted_NS_dividePlots.ipynb` makes plot for paper.
* `ProcrustedEmbedding_<year>.csv` are $n \times d$ dataframes of each year of LAD level data that have been Procrusted as described. 

**19CorrelationMatrices:**
* make correlation matrices of embedding X (use avgX - time averaged X matrix), and factors that I think might be correlated (i.e. driving factors) that I have data for about each nodes.
  
**20PermTestProcrustedData**:
* `permTests_methodA.ipynb` runs permutation tests on the method A way of obtaining the $X_{early}$ and $X_{late}$ Procrusted embeddings. 

**21RobustnessSimulations:**
* `robustnessSimulations.ipynb` contains the simulations given in the Appendix that shows log and Winsorizing make the method more robust (this is that the permutation test is more robust) -i.e. show it works via an example. Consider a student $t$-distribution for heavy tailed data, show difference in the power of the test. Show robustness.
* `robustnessSimulations_taake2.ipynb` - contains the simulations given in the Appendix that shows log and Winsorizing make the method more robust (this is that the permutation test is more robust) -i.e. show it works via an example. Consider a student $t$-distribution for heavy tailed data, show difference in the power of the test. This is the code and results for the paper. Use a KS test to see if transforming the data makes the p-values a closer fit to the Unif[0,1] distribution they should follow.

**22AggregrationJustifications:**
* `AggregationJustifications2.ipynb` -- want to justify that aggregating itself and the choice of level of aggregation are justified. Going to use $z$ scores to assess at the MSOA, LAD, and region level. The findings of this are given in the Appendix to show why we chose LAD level as our choice of data aggregation. 
* `makeA_RegxReg_GB.ipynb` - creates RegxReg matrices from the GB LAD level ones in the folder.

