# mavis3

MAVIS workflow, annotation of structural variants. An application framework for the rapid generation of structural variant consensus, able to visualize the genetic impact and context as well as process both genome and transcriptome data.

## Overview

## Dependencies

* [mavis 3.1.0](http://mavis.bcgsc.ca/)


## Usage

### Cromwell
```
java -jar cromwell.jar run mavis3.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`outputFileNamePrefix`|String|Sample identifier, which will be used for final naming of output files
`reference`|String|The genome reference build. for example: hg19, hg38
`diseaseStatus`|String|Tissue status. For example: diseased
`inputBams`|Array[bamData]|Collection of alignment files with metadata
`svData`|Array[svData]|Collection of SV calls with metadata


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`generateConfig.drawFusionsOnly`|Boolean?|False|flag for MAVIS visualization control
`generateConfig.minClustersPerFile`|Int?|100|Determines the way parallel calculations are organized
`generateConfig.uninformativeFilter`|Boolean?|True|If enabled, only interested in events inside genes, speeds up calculations
`generateConfig.filterMinFlankingReads`|Int?|10|Minimum number of flanking pairs for a call by flanking pairs
`generateConfig.filterMinLinkingSplitReads`|Int?|1|Minimum number of linking split reads for a call by split reads
`generateConfig.filterMinRemappedReads`|Int?|5|Minimum number of remapped reads for a call by contig
`generateConfig.filterMinSpanningReads`|Int?|5|Minimum number of spanning reads for a call by spanning reads
`generateConfig.filterTransHomopolymers`|Boolean?|False|When enabled, transcript sequences containing homopolymer regions are removed
`generateConfig.jobMemory`|Int|6|Memory allocated for this job
`generateConfig.timeout`|Int|6|Timeout in hours, needed to override imposed limits
`runMavis.jobMemory`|Int|36|Memory allocated for this job
`runMavis.timeout`|Int|24|Timeout in hours, needed to override imposed limits


### Outputs

Output | Type | Description
---|---|---
`summary`|File|File with copy number variants, native varscan format
`drawings`|File|Plots generated with MAVIS, collected into a single tar.gz archive
`nscvWT`|File?|Whole transcriptome non-synonymous coding variants. The output file is only generated if variants are found
`nscvWG`|File?|Whole genome non-synonymous coding variants. The output file is only generated if variants are found


## Commands
 This section lists command(s) run by MAVIS 3 workflow
 
 * Running MAVIS 3
 
 #### Generate input configuration file for MAVIS 3
 
 ```
     import json
     import os
 
     ## Set appropriate reference paths
     if REFERENCE == "hg19":
         root =  str(os.environ['HG19_ROOT'])
         mavisRoot = str(os.environ['HG19_MAVIS_ROOT'])
     elif REFERENCE == "hg38":
         root =  str(os.environ['HG38_ROOT'])
         mavisRoot = str(os.environ['HG38_MAVIS_ROOT'])
 
     ## Convert WDL booleans to python booleans
     drawFusionsOnlyPython = eval(DRAW_FUSIONS)
     uninformativeFilterPython = eval(UNINFORMATIVE_FILTER)
     filterTransHomopolymersPython = eval(FILTER_TRANS_HOMOPOLYMERS)
 
     ## Separate input arrays
     b = (sep=' ' BAMS)
     bams = b.split()
     l = (sep=' ' BAM_LIBRARY_DESIGNS)
     bamLibraryDesigns = l.split()
     s = (sep=' ' SV_FILES)
     svFiles = s.split()
     w = (sep=' ' WORKFLOW_NAMES)
     workflowNames = w.split()
     sl = (sep=' ' SV_LIBRARY_DESIGNS)
     svLibraryDesigns = sl.split()
 
     ## Check that appropriate inputs have been supplied for WG and WT analyses
     if ("WG" in bamLibraryDesigns and "WG" in svLibraryDesigns) or ("WT" in bamLibraryDesigns and "WT" in svLibraryDesigns):
         inputs = True
 
     if inputs != True:   
         print("Missing inputs for whole genome and whole transcriptome analyses. Please ensure complete inputs are "
               "supplied for at least one of these analyses.")
 
     else:
         jsonDict = {
             "annotate.draw_fusions_only": DRAW_FUSIONS,
             "cluster.min_clusters_per_file": MIN_CLUSTERS,
             "cluster.uninformative_filter": UNINFORMATIVE_FILTER,
             "summary.filter_min_flanking_reads": FILTER_MIN_FLANKING,
             "summary.filter_min_linking_split_reads": FILTER_MIN_LINKING,
             "summary.filter_min_remapped_reads": FILTER_MIN_REMAPPED,
             "summary.filter_min_spanning_reads": FILTER_MIN_SPANNING,
             "summary.filter_trans_homopolymers": FILTER_TRANS_HOMOPOLYMERS,
             "output_dir": OUTPUT_DIRECTORY,
             "reference.aligner_reference": [
                 ALIGNER_REFERENCE
             ],
             "reference.annotations": [
                 ANNOTATIONS
             ],
             "reference.dgv_annotation": [
                 DGV_ANNOTATIONS
             ],
             "reference.masking": [
                 MASKING
             ],
             "reference.reference_genome": [
                 REFERENCE_GENOME
             ],
             "reference.template_metadata": [
                 TEMPLATE_METADATA
             ]
         }
 
         for index, name in enumerate(workflowNames):
             if name.lower() == "delly":
                 jsonDict["convert"] = {
                     "delly": {
                         "assume_no_untemplated": True,
                         "file_type": "delly",
                         "inputs": [
                             DELLY_FILE_PATH
                         ]
                     }
                 }
             if name.lower() == "starfusion":
                 jsonDict["convert"] = {
                     "starfusion": {
                         "assume_no_untemplated": True,
                         "file_type": "starfusion",
                         "inputs": [
                             STARFUSION_FILE_PATH
                         ]
                     }
                 }
             if name.lower() == "arriba":
                 jsonDict["convert"] = {
                     "arriba": {
                         "assume_no_untemplated": True,
                         "file_type": "arriba",
                         "inputs": [
                             ARRIBA_FILE_PATH
                         ]
                     }
                 }
 
   
         for index, bam in enumerate(bams):
             if bamLibraryDesigns[index] == "WG":
                 jsonDict["libraries"] = {
                     "WG." + OUTPUT_FILE_NAME_PREFIX: {
                         "assign": [
                             "delly"
                         ],
                         "bam_file": BAM_FILE_PATH,
                         "disease_status": DISEASE_STATUS,
                         "protocol": "genome"
                     }
                 }
             if bamLibraryDesigns[index] == "WT":
                 jsonDict["libraries"] = {
                     "WT." + OUTPUT_FILE_NAME_PREFIX: {
                         "assign": [],
                         "bam_file": BAM_FILE_PATH,
                         "disease_status": DISEASE_STATUS,
                         "protocol": "transcriptome",
                         "strand_specific": True
                     }
                 }
 
         for name in workflowNames:
             if name.lower() == "starfusion":
                 jsonDict["libraries"]["WT." + OUTPUT_FILE_NAME_PREFIX]["assign"].append("starfusion")
             if name.lower() == "arriba":
                 jsonDict["libraries"]["WT." + OUTPUT_FILE_NAME_PREFIX]["assign"].append("arriba")
 
         with open("config.json", 'w') as jsonFile:
             json.dump(jsonDict, jsonFile)
 
   ```

 #### Run MAVIS 3

 ```
     
     snakemake --jobs 100 --configfile=CONFIG_FILE -s Snakefile
 
 
     if [ -f OUTPUT_DIRECTORY/summary/MAVIS.COMPLETE ]; then
         ### create an empty zip file, which will be updated with drawings and legends.  if there are none, than the empty file is provisioned out
         echo | zip -q > OUTPUT_FILE_NAME_PREFIX.mavis_drawings.zip && zip -dq OUTPUT_FILE_NAME_PREFIX.mavis_drawings.zip -
 
         ### find all drawing directories, recursively add the drawings
         for draw_dir in `ls -d OUTPUT_DIRECTORY/*OUTPUT_FILE_NAME_PREFIX/annotate/*/drawings`
         do
           zip -qjur OUTPUT_FILE_NAME_PREFIX.mavis_drawings.zip DRAWINGS_DIRECTORY
         done
 
         ### there should be a single mavis_summary_all files
         cp OUTPUT_DIRECTORY/summary/mavis_summary_all_*.tab OUTPUT_FILE_NAME_PREFIX.mavis_summary.tab
 
         ### non-synonymous coding variants are separate into WG or WT files; each may or may not be produced
         if [ -e OUTPUT_DIRECTORY/summary/mavis_summary_WG.*_non-synonymous_coding_variants.tab ];then
           cp OUTPUT_DIRECTORY/summary/mavis_summary_WG.*_non-synonymous_coding_variants.tab OUTPUT_FILE_NAME_PREFIX.WG_non-synonymous_coding_variants.tab
         fi
         if [ -e output_dir_full/summary/mavis_summary_WT.*_non-synonymous_coding_variants.tab ];then
           cp OUTPUT_DIRECTORY/summary/mavis_summary_WT.*_non-synonymous_coding_variants.tab OUTPUT_FILE_NAME_PREFIX.WT_non-synonymous_coding_variants.tab
         fi		  
         exit 0
     fi
     echo "MAVIS job finished but THERE ARE NO RESULTS"
     exit 1
 
   ```
   
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
