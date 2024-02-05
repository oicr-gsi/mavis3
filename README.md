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


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`filterDelly.jobMemory`|Int|24|Memory allocated for this job
`filterDelly.timeout`|Int|6|Timeout in hours, needed to override imposed limits
`generateConfig.drawFusionsOnly`|Boolean?|None|flag for MAVIS visualization control
`generateConfig.minClustersPerFile`|Int?|None|Determines the way parallel calculations are organized
`generateConfig.uninformativeFilter`|Boolean?|None|If enabled, only interested in events inside genes, speeds up calculations
`generateConfig.filterMinFlankingReads`|Int?|None|Minimum number of flanking pairs for a call by flanking pairs
`generateConfig.filterMinLinkingSplitReads`|Int?|None|Minimum number of linking split reads for a call by split reads
`generateConfig.filterMinRemappedReads`|Int?|None|Minimum number of remapped reads for a call by contig
`generateConfig.filterMinSpanningReads`|Int?|None|Minimum number of spanning reads for a call by spanning reads
`generateConfig.filterTransHomopolymers`|Boolean?|None|When enabled, transcript sequences containing homopolymer regions are removed
`generateConfig.jobMemory`|Int|6|Memory allocated for this job
`generateConfig.timeout`|Int|6|Timeout in hours, needed to override imposed limits
`runMavis.jobMemory`|Int|96|Memory allocated for this job
`runMavis.timeout`|Int|72|Timeout in hours, needed to override imposed limits


### Outputs

Output | Type | Description
---|---|---
`summary`|File|File with copy number variants, native varscan format
`drawings`|File|Plots generated with MAVIS, collected into a single tar.gz archive
`nscvWT`|File?|Whole transcriptome non-synonymous coding variants. The output file is only generated if variants are found
`nscvWG`|File?|Whole genome non-synonymous coding variants. The output file is only generated if variants are found


## Commands
 This section lists command(s) run by MAVIS3 workflow
 
 * Running MAVIS3
 
 ```
 
     python3<<CODE 
     import subprocess
     import os
     import shlex
     import gzip
     
     #Separate input arrays
     s = "~{sep=' ' svFiles}"
     svFiles = s.split()
     w = "~{sep=' ' workflowNames}"
     workflowNames = w.split()
     sl = "~{sep=' ' svLibraryDesigns}"
     svLibraryDesigns = sl.split()
 
 
     svFiles_escaped = [shlex.quote(os.path.abspath(path)) for path in svFiles]
 
     for index, name in enumerate(workflowNames):
         if name.lower() == "delly":
             original_delly = svFiles_escaped[index]
             with gzip.open(original_delly, 'r') as f:
               lines = sum(1 for line in f)
 
             if lines > 10000:
                 #Check if other SV callers exist or else survivor can't be run
                 if len(svFiles) > 1:
                     #Run megafusion
                     for index, name in enumerate(workflowNames):
                         if name.lower() == "arriba":
                             arriba_command = f'python3 ~{megafusion_executable} --json ~{megafusion_arriba} --fusion {svFiles_escaped[index]} > arriba.vcf'
                             subprocess.run(arriba_command, shell=True)
                         if name.lower() == "starfusion":
                             starfusion_command = f'python3 ~{megafusion_executable} --json ~{megafusion_starfusion} --fusion {svFiles_escaped[index]} > starfusion.vcf'
                             subprocess.run(starfusion_command, shell=True)
 
                     #Create a copy of the original delly file and increase quality scores to be very high
                     subprocess.run(['cp', original_delly, 'copied_delly.vcf.gz'])
                     subprocess.run(['gunzip', 'copied_delly.vcf.gz'])
                     with open('copied_delly.vcf', 'r') as vcf_file:
                         lines = vcf_file.readlines()
                     with open('copied_delly.vcf', 'w') as vcf_file:
                         for line in lines:
                             if line.startswith('#'):
                                 vcf_file.write(line)
                             else:
                                 fields = line.split('\t')
                                 fields[5] = '1000'
                                 vcf_file.write('\t'.join(fields))
 
                     #Create text file with the modified Delly file and other svFiles
                     input_file_path = 'survivor_input.txt'
                     with open(input_file_path, 'w') as input_file:
                         input_file.write(f'copied_delly.vcf\n')
 
                         #Check for the existence of Arriba and Starfusion files
                         for name in ['arriba', 'starfusion']:
                             sv_file = f'{name}.vcf'
                             if os.path.exists(sv_file):
                                 input_file.write(f'{sv_file}\n')
 
                         #Add any other SV files that are vcfs (e.g. GRIDSS)
                         for sv_file in svFiles:
                             if sv_file.lower().endswith(".vcf"):
                                 input_file.write(f'{sv_file}\n')
 
                     #Run survivor     
                     survivor_command = f'"~{survivor_executable}" merge "survivor_input.txt" 1000 2 0 0 0 1 merged.vcf'
                     result = subprocess.run(survivor_command, shell=True)
                     if result.returncode != 0:
                         raise Exception(f"Error: Survivor command failed with return code {result.returncode}")
 
                     #Look for matching variants from the merged vcf and the original delly file
                     bedtools_command = [
                         'bedtools', 'intersect',
                         '-a', original_delly,
                         '-b', 'merged.vcf',
                         '-header', '-wa', '-u',
                     ]
 
                     subprocess.run(bedtools_command, stdout=open('matched_entries.vcf', 'w'))
 
                 #Filter delly for quality
                 subprocess.run(['bcftools', 'view', '-i', 'FILTER="PASS"', '-O', 'z', original_delly, '-o', 'filtered_delly.vcf.gz'])
 
                 #Add matched variants that were filtered out, to the filtered delly
                 if os.path.exists('merged.vcf'):
                  
                     bedtools_command2 = [
                         'bedtools', 'intersect',
                         '-a', 'matched_entries.vcf',
                         '-b', 'filtered_delly.vcf.gz',
                         '-header', '-v',
                     ]
 
                     subprocess.run(bedtools_command2, stdout=open('variants_not_in_filtered_delly.vcf', 'w'))
 
                     with open('~{outputFileNamePrefix}_mavis_delly.vcf', 'w') as updated_vcf:
                         subprocess.run(['zcat', 'filtered_delly.vcf.gz'], stdout=updated_vcf)
                         subprocess.run(['grep', '-v', '^#', 'variants_not_in_filtered_delly.vcf'], stdout=updated_vcf)
 
                 else:
                     subprocess.Popen(['gunzip', '-c', 'filtered_delly.vcf.gz'], stdout=open('~{outputFileNamePrefix}_mavis_delly.vcf', 'w'))
 
             else:
                 subprocess.Popen(['gunzip', '-c', svFiles_escaped[index]], stdout=open('~{outputFileNamePrefix}_mavis_delly.vcf', 'w'))
 
     CODE
   ```
 ```
 
     ## Use python snippet to generate config file
     python3<<CODE 
     import json
     import os
 
     #Set appropriate reference paths
     if "~{reference}" == "hg19":
         root =  str(os.environ['HG19_ROOT'])
         mavisRoot = str(os.environ['HG19_MAVIS_ROOT'])
     elif "~{reference}" == "hg38":
         root =  str(os.environ['HG38_ROOT'])
         mavisRoot = str(os.environ['HG38_MAVIS_ROOT'])
 
     #Convert WDL booleans to python booleans
     drawFusionsOnlyPython = eval("~{drawFusionsActual}".title())
     uninformativeFilterPython = eval("~{uninformativeFilterActual}".title())
     filterTransHomopolymersPython = eval("~{filterTransActual}".title())
 
     #Separate input arrays
     b = "~{sep=' ' bams}"
     bams = b.split()
     l = "~{sep=' ' bamLibraryDesigns}"
     bamLibraryDesigns = l.split()
     s = "~{sep=' ' svFiles}"
     svFiles = s.split()
     w = "~{sep=' ' workflowNames}"
     workflowNames = w.split()
     sl = "~{sep=' ' svLibraryDesigns}"
     svLibraryDesigns = sl.split()
 
     #Check that appropriate inputs have been supplied for WG and WT analyses
     if ("WG" in bamLibraryDesigns and "WG" in svLibraryDesigns) or ("WT" in bamLibraryDesigns and "WT" in svLibraryDesigns):
         inputs = True
 
     if inputs != True:   
         print("Missing inputs for whole genome and whole transcriptome analyses. Please ensure complete inputs are "
               "supplied for at least one of these analyses.")
 
     else:
         jsonDict = {
             "annotate.draw_fusions_only": drawFusionsOnlyPython,
             "cluster.min_clusters_per_file": ~{minClustersActual},
             "cluster.uninformative_filter": uninformativeFilterPython,
             "summary.filter_min_flanking_reads": ~{filterMinFlankingActual},
             "summary.filter_min_linking_split_reads": ~{filterMinLinkingActual},
             "summary.filter_min_remapped_reads": ~{filterMinRemappedActual},
             "summary.filter_min_spanning_reads": ~{filterMinSpanningActual},
             "summary.filter_trans_homopolymers": filterTransHomopolymersPython,
             "output_dir": "output_dir_full",
             "reference.aligner_reference": [
                 mavisRoot+"~{alignerReference}"
             ],
             "reference.annotations": [
                 mavisRoot+"~{annotations}"
             ],
             "reference.dgv_annotation": [
                 mavisRoot+"~{dgvAnnotation}"
             ],
             "reference.masking": [
                 mavisRoot+"~{masking}"
             ],
             "reference.reference_genome": [
                 root+"~{referenceGenome}"
             ],
             "reference.template_metadata": [
                 mavisRoot+"~{templateMetadata}"
             ],
             "convert": {},
             "libraries": {}
         }
 
         for index, name in enumerate(workflowNames):
             if name.lower() == "delly":
                 if "delly" not in jsonDict["convert"]:
                     jsonDict["convert"]["delly"] = {
                         "assume_no_untemplated": True,
                         "file_type": "delly",
                         "inputs": [
                             str(svFiles[index])
                         ]
                     }
 
             if name.lower() == "starfusion":
                 if "starfusion" not in jsonDict["convert"]:
                     jsonDict["convert"]["starfusion"] = {
                         "assume_no_untemplated": True,
                         "file_type": "starfusion",
                         "inputs": [
                             str(svFiles[index])
                         ]
                     }
 
             if name.lower() == "arriba":
                 if "arriba" not in jsonDict["convert"]:
                     jsonDict["convert"]["arriba"] = {
                         "assume_no_untemplated": True,
                         "file_type": "arriba",
                         "inputs": [
                             str(svFiles[index])
                         ]
                     }
   
         for index, bam in enumerate(bams):
             if bamLibraryDesigns[index] == "WG":
                 if "WG." + "~{outputFileNamePrefix}" not in jsonDict["libraries"]: 
                     jsonDict["libraries"]["WG." + "~{outputFileNamePrefix}"] = {
                         "assign": [],
                         "bam_file": bams[index],
                         "disease_status": "~{diseaseStatus}",
                         "protocol": "genome"
                     }
             if bamLibraryDesigns[index] == "WT":
                 if "WT." + "~{outputFileNamePrefix}" not in jsonDict["libraries"]:
                     jsonDict["libraries"]["WT." + "~{outputFileNamePrefix}"] = {
                         "assign": [],
                         "bam_file": bams[index],
                         "disease_status": "~{diseaseStatus}",
                         "protocol": "transcriptome",
                         "strand_specific": True
                     }
 
         for name in workflowNames:
             if name.lower() == "delly":
                 jsonDict["libraries"]["WG." + "~{outputFileNamePrefix}"]["assign"].append("delly")        
             if name.lower() == "starfusion":
                 jsonDict["libraries"]["WT." + "~{outputFileNamePrefix}"]["assign"].append("starfusion")
             if name.lower() == "arriba":
                 jsonDict["libraries"]["WT." + "~{outputFileNamePrefix}"]["assign"].append("arriba")
 
         with open("config.json", 'w') as jsonFile:
             json.dump(jsonDict, jsonFile)
 
     CODE
   ```
 ```
 
     set -eu -o pipefail
     
     snakemake --jobs 40 --configfile=~{configFile} -s /.mounts/labs/gsiprojects/gsi/Investigations/GRD-603.Mavis3/Delly_filtering/Cleaning_Delly_filtering/3Attempt/OCT_011303/Updated_Snakefile/Snakefile
 
 
     if [ -f output_dir_full/summary/MAVIS.COMPLETE ]; then
         ### create an empty zip file, which will be updated with drawings and legends.  if there are none, than the empty file is provisioned out
         echo | zip -q > ~{outputFileNamePrefix}.mavis_drawings.zip && zip -dq ~{outputFileNamePrefix}.mavis_drawings.zip -
 
         ### find all drawing directories, recursively add the drawings
         for draw_dir in `ls -d output_dir_full/*~{outputFileNamePrefix}/annotate/*/drawings`
         do
           zip -qjur ~{outputFileNamePrefix}.mavis_drawings.zip $draw_dir
         done
 
         ### there should be a single mavis_summary_all files
         cp output_dir_full/summary/mavis_summary_all_*.tab ~{outputFileNamePrefix}.mavis_summary.tab
 
         ### non-synonymous coding variants are separate into WG or WT files; each may or may not be produced
         if [ -e output_dir_full/summary/mavis_summary_WG.*_non-synonymous_coding_variants.tab ];then
           cp output_dir_full/summary/mavis_summary_WG.*_non-synonymous_coding_variants.tab ~{outputFileNamePrefix}.WG_non-synonymous_coding_variants.tab
         fi
         if [ -e output_dir_full/summary/mavis_summary_WT.*_non-synonymous_coding_variants.tab ];then
           cp output_dir_full/summary/mavis_summary_WT.*_non-synonymous_coding_variants.tab ~{outputFileNamePrefix}.WT_non-synonymous_coding_variants.tab
         fi		  
         exit 0
     fi
     echo "MAVIS job finished but THERE ARE NO RESULTS"
     exit 1
 
   ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
