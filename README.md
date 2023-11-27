# GB-Hyperlinks_N-S-Divide
Data and Code Repository to accompany the paper: 'Assessing the North-South Divide in Britain’s Digital Landscape through Stable Dynamic Embedding of Spatial Web Data'. 

The datasets used are:
*the Geoindex and Host Link Graphs which can be downloaded from JISC UK Web Domain Dataset: `https://data.webarchive.org.uk/opendata/ukwa.ds.2/`. This repository contains the code used to process the files, but not the files themselves. 
* ONS postcode lookup, which has many versions that can be downloaded from: https://geoportal.statistics.gov.uk/search?q=postcode%20lookup. (This project accessed a Nov2021 version)


**01Bash:**
* ```MSOA_REG11NM_lookup_bash.sh``` outputs a file called ```MSOA_REG11NM_lookup.csv``` which has a column for MSOA and Region. It uses the file `Output_Area_to_LSOA_to_MSOA_to_Local_Authority_District_(December_2017)_Lookup_with_Area_Classifications_in_Great_Britain.csv` downloaded from: https://geoportal.statistics.gov.uk/datasets/national-statistics-postcode-lookup-november-2021/about
* The same file from the ONS is used to make a postcode (pcds) to MSOA (msoa11) lookup, and is done via the script `MSOA_all_nov2021_bash.sh` and the output is called `MSOA_UK_nov2021.csv`.
* `lookup_bash.sh` takes in `all.1.csv` to get out the host and pc only called `lookup_all.csv` -- `all.1.csv` is a file from Emmanouil's previous project of all archived websites with one postcode. 

**02DataProcessing**: Script to make lookup table (host/postcode)
* `lookup_all_script.R` uses output `host-pc.csv` to make lookup table called `lookup_all.csv`. Keeps only one unique pc for a host.

**04HPCScripts:** 
* `MSOA_1996-2003_2005-2010_filter_couk_couk.sh` runs the R script called `1996-2003_2005-2010_couk_couk_filter.R` to filter for each year in the linkage data, links between 2 ".co.uk" hosts. THe output files take the form `<year>-couk-couk-linkage.tsv`.
* * ```MSOA_2005-2010_make_adjacencyMatrices_v2.sh``` uses the R script ```MSOA_2005-2010_MakeAdjMatr_v2.R``` , which has 2 different filters to remove self-hyperlinks (within the website). It does this by comparing host1 and host2. It outputs A_<year>_v2.mtx files for 2005-2010 inclusive. Changes made the filters in this document on 22/04/2022. This takes the `<year>-couk-couk-linkage.tsv` files and the MSOA-postcode lookup and matches both sides of the *'.co.uk'* to *'.co.uk'* links to an MSOA via their postcode, then aggregates this data by MSOA to output matrices `A_for_<year>.mtx` files, where $A_{ij} represents the number of hyperlinks from MSOA[i] to MSOA[j] -- so these are weighted, directed matrices. Running script also prints an output of general summary of the A matrices to the console. 

**05AdjMatrices:**
* contains the A (adjacency matrices) for 1996-2010.
* Script ```EDA_on_A_matrices.R``` which goes through the A matrices (as A(i,j)+A(j,i) to symmetrise), and sees some things about the connected components.
* frequency tables, and largest connected component csv files - tell you information about these graphs. 2004 data so large it is omitted in most work in this section. Can find nice summaries. File includes plots of the LCC in each year. 
* Can see histograms of how many of the 14 LCCs for 1996-2003_2005-2010 an MSOA appears in. Also made a histogram for 2005-2010.
* The LCC for 2005-2010 contains more than half the nodes in the graph for all 6 of these years. 
* ```MSOA_node_table.R``` script output a csv file called *nodes_MSOA.csv* of MSOA and index (which maps an MSOA to the index of the A matrices).

**06Maps**:
* ```EDA_A_outgoing_2005-2010_map.Rmd``` takes in the A matrices, and formats into a vector of "how many outgoing hyperlinks from each MSOA". Then this is used with a MSOA ```geojson``` shapefile, to plot these values on a map (images of such are in folder). Then want to put into a hexagon shaped cartogram, to try and display this data in an easily interpretable way.
* ```check_hyperlinks_values.Rmd```: could see that some values in A matrices were very large. This investigates this, and highlights the majority of entries in A that are large are within an MSOA. This leads to a **v2** A matrices script in the **04HPCScripts** folder, that removes any data entries between a host and themselves. I will see if this is sufficient and reasonable. 20/04/2022 - also want to consider a self-defined cap on the largest value for number of hyperlinks in the data, before we exclude the data entry, as we think it is a directory like site in the data. 
* ```EDA_A_London_2005-2010.Rmd```. Want to subset the A matrices data into just MSOAs that are in *London*, and make similar plots as for the whole of E&W (map by MSOA and a hex map), as well as a plot that shows where the links are going to for a given MSOA. 
  
  
**07SpectralEmbeddingsOfA**:
* ```UASE_on_A_v2_2005-2010.Rmd``` has initial partial SVD of unfolded A, and a scree plot to go with.
Then move on to use Python via Jupyter lab!!
* ```uase_with_ed.ipynb``` does a spectral embedding, and degree correction. It makes plots that are coloured by
  1. is or is not in London for all data
  2. for England only, which LEP is the MSOA in 
* ```UASE_on_A_v2_2005-2010_NUTS.ipynb``` - trying to find a way to colour data for all MSOAs in the dataset.
  1. By Rural-Urban status (8 factor levels) only for E&W
  2. By region for all of GB
  
**08LAD**:
* `LAD_level_investigation.ipynb` is a notebook to explore the outwards A matrices that are at a LAD x LAD level. See whether embeddings work at this granularity. Also a bit of incomplete workings to do with embedding the Reg x LAD matrices.
* `DanCode_plusExtra.R` is some code that uses a dendrogram and hierarchical clustering to order the LADs for the heatmap plots. If you add the 6 A matrices for all years together to make a "mean" matrix and take the transform $(A + A^T)/2$ for the "mean" matrix, you can see small clusters appear along the diagonal. It also outputs some plots from running an svd on the matrices and the "mean" matrix, which show geographical structure in some ways. 
* Folder contains animations of the LAD x LAD embeddings, and lots of the heatmaps which are appropriately labelled. 
* The heatmaps for each year highlight some LADs in 2005, 2006, and 2008 seem overtly highly connected, so investigating these is in **09**. 
  
**09OddLADInvestigation**:
* `Odd_LAD_explore.ipynb` and `oddLADexplore.R` - using both at the same time to try and get to the root of causes of the highly outwardly connected LADs seen in some years. Start with 2006, where Wolverhampton and Leicester are highly connected. 
* `potential_odd_websites2006.csv` contains all the websites and the Postcodes that could be weird, so that this can be joined with the `2006-couk-couk-linkage.tsv` file to find out. 
* `potential_odd_hyperlinks2006_HPC.R` and the `.sh` file are to run the join in HPC and output a csv file of all the hyperlink data entries that are associated with host1 being in the Wolverhampton or Leicester LAD. 
* `<year>_topXXOddHyperlinks.csv` comes from `oddHyperlinkRxplore<year>.R` script, and looks at the top XX odd websites. The csv files are to use with `oddLADexplore.R where it is noted which postcodes are odd.
* `oddHyperlinks_whereTo.R` is a script to find all the raw data lines (from the hyperlinks data stored in the **<year>-couk-couk** files), where the host1 is potentially odd, to allow visual examination of the host2 they are linking to, to understand if they are legitimate, and what is going on. Uses the **potential_odd_hyperlinks<year>.csv** files and filters these. 
  
**10LADDegreeDistributions:**
* `histograms2005-2010.R` just experiment plotting the values of the **number of hyperlinks** between sites on a histogram. The point is to see at what maximum values cutting off the top XX% helps. 
* The *1PC* files are more useful, as this just looks at the observations in the raw data that go towards to A matrices (which are for 1PC only).
  
**11LADEmbeddingsCroppedData:**
* Makes A_LAD matrices, that have had the top 0.01% of highest number of hyperlinks observations (per year) removed from the raw data, before they are aggregated to adjacency matrices. This is **version3**.
  
**12BinaryHyperlinkConnections:**
* Here, instead of caring and aggregating **how many hyperlinks between places**, we are making matrices of **how many website to website connections between places** are there. This means replacing the *# hyperlinks* by a 1 if non-zero and a 0 if zero, to make an adjacency matrix. Can think of this as taking the huge adj. matrix of website level and replacing all non-zero values with a 1. Then they are aggregated to MSOA and then LAD level by summing. `MSOA_2005-2010_makeAdjMatr_binary.R` does this, then there is a Python script `makeA_LAD_binary.ipynb`, that aggregates the MSOA level matrix to LAD level. (Also Winsorizing the matrices and log10 on them is done). This makes **version 4 of the MSOA x MSOA matrices**.
* Histograms of various degree distributions, from log10 scale, to adding also Winsorization at (0, 0.99) level. Heatmaps of the degree distributions also. These help to indentify appropriate transformations of the data to make less heavy tailed.
* There is then `2D_colourmap_reference.pdf` which shows a 2D rgb scale that E&W LADs are then matched to by geography: `FilledMapLADS_E_W.pdf`.
* Assortment of SVD (embedding) plots then, where titles suggest transformations and the points are labelled by geography. Also a scree plot to see appropriate number of dimensions to consider.
* `Animated_plots_SVDbinaryLog10Wins1.R` plots using the title name transformations, and also creates the `animationData.csv` file which is used again.
* `animatedData_makingIt.ipynb` takes in the csv `animationData.csv` and creates animations across all years of the points in PC1/PC2 and PC3/PC4, where points are coloured by geographical location. There is also an animation of PC1/PC2 for 2009-2010 which is coloured by **number of active enterprises**.
* Uses `Annual_pay-Gross-all_employee_jobsaUK2010.csv` from https://www.ons.gov.uk/employmentandlabourmarket/peopleinwork/earningsandworkinghours/datasets/placeofresidencebylocalauthorityashetable8 to plot the PCs by pay in `animationData_makingIt.ipynb` notebook.
* `transformationsOfAPlot.R` generates plots to show what we do to the E&W LAD level A matrices at this stage (See **16** for updated plots for GB where we conc then Wins).
**12-1-histogram_weights:**
* histograms of the weights in the non-binary A matrices. These are to justify using the binary A to reduce the heavy tails of the distributions of A.
**13DistanceMeasures:**
* `create_mean_in_degree_csv.R` and an `out` file too create 2 csv, where they use all the `A_LAD_<year>_binary.mtx` files to find out the mean inward degree for each LAD. This is to be used to order the heatmaps.
* `distanceMeasuring.ipynb` looks at the cosine distance of points (LADs) in the embeddings, and looks at how much these change between consecutive years. Heatmaps of these are made, where we order to try and understand which points are moving a lot and which are not - we intend to find some understanding of *why* nodes move a lot or do not. 
* `Animated_plots_SVDbinaryLog10Wins1.R` is the file that creates the `animationData.csv` file. 
* `cosineDistancesMovedLAD.csv` is a file output from the notebook above, which contains the cosine distances moved for each LAD. It is to be used by the `distanceMovedAnalysis.R` file, which reads in this `csv` file and tries to cluster them to find out any patterns about how the points move. - hclust not so useful.
**14NorthSouthDivide:**
* **LADS_list.csv** tells you the order of the indexes of the LAD x LAD matrices
* `Regions_list.csv` tells you the order of the indexes of the Reg x Reg matrices
* `NorthSouthDivide.ipynb` converts the LAD x LAD matrices to Reg x Reg matrices and saves them. It then goes on to create embeddings at the LAD level, which cannot be directly compared to the LAD level UASE embeddings, but may help to provide insight into the behaviour N/S.
    End of the file there is some code colouring the LAD level embeddings by urban-rural classification, using a file to map the categories found in 08 ( copied in).
* `permutation_test_N-S-divide_short.ipynb` runs permutation tests, as has visual outputs for this. It also looks at the year-on-year embeddings average position for the N/S to detect whether there is trend. There are graphs from this that provide justification for us permutation testing the hypothesis that we do. 
* `permuation_testing_other_things.ipynb` this tests other hypotheses using the same methods as the N/S divide test. Lots of the code is not relabelled all the way through, but is labelled correctly for the plots and at the start with titles so you can follow what is being followed. 
* `simulations.ipynb` runs a model with a SBM of 2 communities. It shows how the permutation tests look when you know what the outcome should look like. This helps to justify the plots, add explainability when we apply these methods to real data (where we do not know the underlying distribution that the nodes come from). 
  
**15 UKincScot:**
Note: I've called it UK everywhere, but it's actually GB as NI is not included. 

v5 is binary still, just includes GB. 
* `MSOA_LAD_UK_lookup_bash.sh` doesn't actually work. Gets stuck because strings inside cells have commas in them, but can make the lookup file a usable size in python (lucky for me!). -- *so don't try to use*. Outputs file `MSOA_LAD_UK_lookup.csv` to use. 
* `makeA_LAD_v5_UKincScot.ipynb` takes in the file `PCD_OA_LSOA_MSOA_LAD_FEB20_UK_LU.csv` which was downloaded from a link given in the notebook. The start of the notebook makes this file into a MSOA to LAD lookup table for the UK, which it saves as the file `MSOA_LAD_UK_lookup.csv` which has been saved and then can just be used (this is where it is made). The notebook makes the LADxLAD matrices (v5) for the UK, so Eng, Wales + Scot, using the v4 binary MSOA x MSOA matrices made in **12**.
* `LADS_list_UK.csv` is a list of all the LADs in E, W & S. `A_LAD_<year>_UK_v5.mtx` are the matrices for the whole UK at LAD x LAD level.
* `matrixSVD_plots_UK_binaryLog10Wins.R` uses the shp file from: 
https://geoportal.statistics.gov.uk/datasets/ons::local-authority-districts-december-2016-generalised-clipped-boundaries-in-the-uk-1/explore?location=53.588449%2C-0.820715%2C5.52
and does the SVD on the v5 A LADxLAD matrices. Makes the "animationData" and saves this as `animationDataUK.csv`. 
Second part of the `R` script makes plots: scree plots, 2D reference map plots, embedding plots. 
* `permutationTestsUK.ipynb` script to run permuation tests. (Also has code to output a `LAD_NorthIndicator_GB_lookup.csv` file, which is copied into *18* and used there)
* `animations.ipynb` makes some animation data, does some plots of the count of active enterprises data, and outputs the CLEAN version of the active ent data (which is copied to 18 and used in the 1dim Procrustes stuff).
  
**15-1 WinsorizationPoint:**
* `Winsorizing_beforeAfter_concatenation.R` implements a distance measure for the embeddings and compares for concatenating then Winsorizing and vice versa.
  
**16GB_ConcThenWins:**
* `matrixSVD_plots_UK_binaryLog10ConcWins.R` file takes in the LAD level A matrices, then concatenates, logs, then Winsorizes. Then does SVD to calculate embedding, and makes some plots of the embedding and a scree plot, and the file `animationDataGB_ConcWins.csv`.
* `permutationTestsGB_ConcWins.ipynb` uses the animation data to do the permutation tests for N/S, U/R, L/nL and plots graphs and outputs the values of the p-values. 
* `transformationsAPlots_GB_ConcWins.R` generates similar plots to the **12** file, but does for when we conc then Wins the data, and for GB not just E&W.
  
**17NorthernIrelandData:**
* `PCD_SOA_NI_lookup_bash.sh` uses the file `ONSPD_AUG_2011_UK_O.csv` to crop it down to the PCD and SOA for Northern Ireland, and outputs this as the file `postcode_SOA_NI_lookup.csv`.
* SOA for NI is the equivalent for MSOA, so build up a new R script for the HPC to build binary A matrices at the MSOA x MSOA level: `MSOA_UK_2005-2010_makeAdjMatr_binary.R`. There is a bash script to go with this. They are in the *binaryA folder on HPC as existing files were there and additional files needed for NI have been uploaded. Outputs **A_for_<year>_binary_MSOA_UK.mtx**. These will then need to be fed into python script to make LAD x LAD (LGD is the NI equvialent) matrices. 
* `make_nodes_MSOA_UK_csv_file.R` outputs the file `nodes_MSOA_UK.csv` to use in python script. 
* `makeA_LAD_v6_UK.ipynb` to convert MSOA level A to LAD level A.  
* https://www.nisra.gov.uk/publications/11-dc-look-tables-and-guidance-documents table to go SOA to LGD in NI.
*  `SVD_plots_UK_ConcWins.R` makes animation Data, SVD plots, scree plots, maps filled and outline of the UK.
* `permutationTestsUK_ConcWins.ipynb` uses the `animationDataUK_ConcWins.csv` data to do the permutation tests. 
* In NI use data from https://www.ninis2.nisra.gov.uk/public/PivotGrid.aspx?ds=3635&lh=73&yn=2006,2009&sk=10&sn=People%20and%20Places&yearfilter= to classify the LGD2014 (LAD equivalents as U/R), based on the % of Urban dwellings in the LAD. Saved as `urbanDwellingsPercentage_LGD2014_NI.csv`.
* csv file of *all* permutation tests for different areas and data pre-processing are in the file: `permutationTestPvaluesALL.csv`.
* `howManyUniqueHostsInHostLinkData2005-2010.R` and `avgNumberWebsitesLinkedToInHostLinkData2005-2010.R` is to find out how many websites are in the subset of the host-link graph data we use for the paper. 
* `makeA_RegxLAD_UK.ipynb` makes Region x LAD matrices for the UK: **A_RegxLAD_<year>_UK.mtx**.
* `embedA_regxLAD.ipynb` does an SVD and creates an embedding of these matrices. Tries to output some animations (html) files.

**18Procrustes:**
* `make_A_LAD_GB.ipynb` crops A_v6_LAD to make them just for the GB.
* `ProcrustesGB.R` tries to do a Procrustes transformation to the SVD of A (which is given by X), to rotate it and project it into a lower dimension. It uses an iterative algorithm to do this, to try and fit $\tilde{X}$ to G (geography - lat/long). Then uses the $\tilde{X} \in \mathbb{R}^{n \times 2}$ to give the nodes (LADs) a colouring on a 2D heatmap scale, that are based on their $\tilde{X}$ value, and an opacity based on the "connectivity" of the node. Hope that this looks like geography. 
* `ProcrustesGB_avg.R` is a big improvement on the previous attempt.
* `ProcrustesGB_X_avg_v2.R` has lots of work in it, makes the plots and does the iterative Procrustes transformation. Does the zoom in on London, as well as doing some X_early and X_late Procrustes transformations. - Think I need to extract this information to then be able to apply the permutation tests on the data. 
* `Procust_early_late_for_PermTests.R` is a script to procrust the embedding into $C = (long, lat,  no.enterprises, zeros)$. Method (a) procrustes $X_{all average}$ into $C$ to learn the transformation parameters $s, R, tt$. It then uses these values to Procrust $X_{early}$ and $X_{late}$ into $C$. Method (b) does: $X_{early} \in \mathbb{R}^{n \times 4}$ and $X_{late}$ into $C = (long, lat,  no.enterprises, zeros)$. Then export the transformed dataframes as csv files to do some permutation tests on. (24/04/2023) The output csv files are copied into **20** to run permutation tests on them.
* `CheckNumberOfLinksPerYear_12Sep2023.ipynb` - plot to show that the number of links is the pattern we see in dim1 of the UASE embedding of the data (Figure 5 in the paper)
* `Procrust_NS_divide_plots.R` is done (oct2023) to make the plots for the paper wrt the Procrusted Dimensions rather than the embedding dimensions for better reading. `YhatProcrusted_allYears.csv` contains the Procrusted embeddings for each year all stacked into one matrix. 
* `Procrusted_NS_dividePlots.ipynb` make plot for paper. 

**19CorrelationMatrices:**
* make correlation matrices of embedding X (use avgX - time averaged X matrix), and factors that I think might be correlated (i.e. driving factors) that I have data for about each nodes.
  
**20PermTestProcrustedData**:
* `permTests_methodA.ipynb` runs permutation tests on the method A way of obtaining the $X_{early}$ and $X_{late}$ Procrusted embeddings. 

**21RobustnessSimulations:**
* `robustnessSimulations.ipynb` Need some appendix simulation that shows log and Winsorizing make the method more robust (I think this is that the permutation test is more robust) -i.e. show it works via an example. Consider a student $t$-distribution for heavy tailed data, show difference in the power of the test. Show robustness.
* `robustnessSimulations_taake2.ipynb` - these have new code rather than trying to reuse the old code. This is the code and results for the paper. Use a KS test to see if transforming the data makes the p-values a closer fit to the Unif[0,1] distribution they should follow.

**22AggregrationJustifications:**
* Want to justify that aggregating itself and the choice of level of aggregation are justified. Probably won't do at the website level for now, just becuase it is so computationally expensive. Going to use $Z$ scores to assess at the MSOA, LAD, and region level.
* `makeA_RegxReg_GB.ipynb` - can't find proper RegxReg matrices so will make from the GB LAD level ones in the folder.

