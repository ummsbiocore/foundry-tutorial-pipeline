**RNA-sequencing Pipeline**
Version: 1.2.0




##### **DE module**: Performs differential expression analysis using DESeq2 software (see the [<span style="color:blue">documentation</span>](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html),  the [<span style="color:blue">code</span>](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), and the [<span style="color:blue">paper</span>](https://doi.org/10.1186/s13059-014-0550-8)) or Limma Voom software (see the [<span style="color:blue">documentation</span>](https://bioconductor.org/packages/release/bioc/html/limma.html) and the [<span style="color:blue">paper</span>](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29)) . 

**Replicates are required** for differential expression analysis. At least 3 for each condition is recommended. See [<span style="color:blue">Schurch et.al.</span>](https://dx.doi.org/10.1261%2Frna.053959.115) for more information on how many replicates are necessary.
<br>


##### **Module Input:**
The minimum input requirements are:

1) A "counts" file which contains gene or transcript counts. **Requirements:**
  *  Tab separated (TSV) with a header containing sample names.
  *  The first column must contain features (genes or transcripts).
  *  No hyphens/dashes/spaces are allowed in column names.
  *  Values in the first column (features) must be unique.

2)  A "groups" file which contains sample information. **Requirements:** 
  * Tab separated (TSV) with a header
  * **sample_name** and **group** columns are required.
  * For the sample_name column, please use the names entered while creating a collection.
  * Each sample_name must appear as a column in the header from the counts file.
  * No hyphens/dashes/spaces are allowed in column names.
  * No hyphens/dashes/spaces are allowed in the group column.
  * Additional columns can be included which contain metadata about each sample

3) A "comparison" file which specifies which groups to compare in DE analysis. **Requirements**:  
  * Tab separated (TSV) with a header
  * column names are exactly:
  
   | **controls** | **treats** | **names** | 
	 
  * No hyphens/dashes/spaces are allowed in either the column names or the values.
  * The values in treats and controls columns must exist in the group column from the groups file.
  * The names column will control the prefix of output files.
  * An optional fourth column is allowed called **grouping_column**. This will override the column from the group file that is used for comparisons. This value will also override the design formula such that it becomes "~*grouping_column*". More complicated design formulas are not supported when using the optional fourth column.
  
4) Experimental design formula in R formula notation.
  * This is used in building the statistical model for DE analysis.
  * By default this is set to "~group" which means we are interested in comparing    
    different groups, as listed in the group file. 
  * For most users the default value would be enough and should not be changed.
  * Advanced users can change the default value to other formulae using the column names from the groups file.

Batch removal from PCA and Heatmap plots: Batches are handled internally for differential gene expression by
DESeq using an appropriate design formula (see the examples below). However, this process does not involve batch removal
 from the data. It may be desirable to remove batch effects from the data when creating PCA plots. To do so, enable
 batch correction and set the "Correction column" to a column from the groups file that should be used for batch correciton.
 The Biological variable column should be set to "group" in most cases (when the DESeq design is set to ~group)

#### Example
1) You want to compare two conditions with three replicates each: compound treatment and placebo control

**Groups file:**

| **sample_name** | **group**   |
|--------------------|-------------|  
| sample1                   | placebo        |
| sample2                   | placebo       |
| sample3                   | placebo       |
| sample4                   | compound  |
| sample5                   | compound  |
| sample6                   | compound  |

<br>

**Comparison file:**
 
|** controls** | **treats**   | **names**                            |
|--------------|-------------|----------------------------|
|placebo          | compound   | compound_vs_placebo  |

<br>

#### DESeq2 steps
*  Create a DESeq Data Set (DDS) from input gene counts, filter using the provided minimum sample and read counts values.
*  Output normalized gene expression and z-scores from the DDS for all genes/samples making through the filter (AllDetectedGenes)
*  Create a PCA plot for the first 2 principal components using a variance stabilization transformed version of the gene counts.
*  Perform hierarchical clustering for the specified number of clusters on z-scores and normalized expression values for AllDetectedGenes
*  Create heatmap of the hierarchical clustering results.
*  Perform differential expression analysis for each comparison listed in the compare file
*  For each comparison create volcano and MA plots for DE genes.


#### Adding Sequences To Reference Genome 

1. Turn on the "run_Download_Genomic_Sources" option.
2. Activate the "add_sequences_to_reference" option and configure its settings. You can input your custom sequence in FASTA format. If you prefer to provide GTF file, it is also supported. If you don't provide one, it will be automatically generated.
3. Enable the relevant aligner-specific index-building options (e.g., "build_STAR_index," "build_RSEM_index") based on your chosen aligners. This step ensures the inclusion of custom FASTA/GTF files into the genome.
