version 1.0

struct bamData {
  File bam
  String libraryDesign
}

struct svData {
  File svFile
  String libraryDesign
  String workflowName
}

struct mavisResources {
    String alignerReference
    String annotations
    String dgvAnnotation
    String masking
    String referenceGenome
    String templateMetadata
}


workflow mavis3 {
  input {
    String outputFileNamePrefix
    String reference
    String diseaseStatus
    Array[bamData] inputBams
    Array[svData] svData
  }

  parameter_meta {
    outputFileNamePrefix: "Sample identifier, which will be used for final naming of output files"
    reference: "The genome reference build. for example: hg19, hg38"
    diseaseStatus: "Tissue status. For example: diseased"
    inputBams: "Collection of alignment files with metadata"
    svData: "Collection of SV calls with metadata"
  }

  String sanitizedSID = sub(outputFileNamePrefix, "_", ".")

  ## Setting modules and reference paths
  String filterdelly_modules = "megafusion/0.0.1 survivor/1.0.7 bcftools/1.9 bedtools/2.27" 
  String megafusion_executable = "$MEGAFUSION_ROOT/bin/MegaFusion.py"
  String megafusion_arriba = "$MEGAFUSION_ROOT/arriba.json"
  String megafusion_starfusion = "$MEGAFUSION_ROOT/starfusion.json"
  String survivor_executable = "$SURVIVOR_ROOT/Debug/SURVIVOR"
  
  Map[String,String] mavis_modules_by_genome = { "hg19": "mavis/3.1.0 hg19-mavis/3.1.0 hg19/p13", "hg38" : "mavis/3.1.0 hg38-mavis/3.1.0 hg38/p12" }
  String mavis_modules = mavis_modules_by_genome [ reference ]

  Map[String,mavisResources] resources = {
    "hg19": {
      "alignerReference": "/hg19.2bit",
      "annotations": "/ensembl69_hg19_annotations_with_ncrna.v3.json",
      "dgvAnnotation": "/dgv_hg19_variants.tab",
      "masking": "/hg19_masking.tab",
      "referenceGenome": "/hg19_random.fa",
      "templateMetadata": "/cytoBand.txt"
    },
    "hg38": {
      "alignerReference": "/hg38.2bit",
      "annotations": "/ensembl79_hg38_annotations.v3.json",
      "dgvAnnotation": "/dgv_hg38_variants.tab",
      "masking": "/hg38_masking.tab",
      "referenceGenome": "/hg38_random.fa",
      "templateMetadata": "/cytoBand.txt"
    }
  }


  scatter(b in inputBams) {
    String bams = "~{b.bam}"
    String bamLibraryDesigns = b.libraryDesign
  }

  scatter(s in svData) {
    String tempSvFiles = "~{s.svFile}"
    String workflowNames = s.workflowName
    String svLibraryDesigns = s.libraryDesign
  }

  scatter(s in svData) {
    if (s.workflowName == "delly") {
      call filterDelly {
        input:
          outputFileNamePrefix = sanitizedSID,
          svFiles = tempSvFiles,
          workflowNames = workflowNames,
          svLibraryDesigns = svLibraryDesigns,
          modules = filterdelly_modules,
          megafusionExecutable = megafusion_executable,
          megafusionArriba = megafusion_arriba,
          megafusionStarfusion = megafusion_starfusion,
          survivorExecutable = survivor_executable
      }
    }    
    String svFiles = select_first([filterDelly.mavisDelly,"~{s.svFile}"])
  }

  call generateConfig { 
    input: 
      outputFileNamePrefix=sanitizedSID, 
      diseaseStatus=diseaseStatus, 
      alignerReference=resources[reference].alignerReference, 
      annotations=resources[reference].annotations, 
      dgvAnnotation=resources[reference].dgvAnnotation, 
      masking=resources[reference].masking, 
      referenceGenome=resources[reference].referenceGenome, 
      templateMetadata=resources[reference].templateMetadata,
      reference=reference,
      modules=mavis_modules,
      bams=bams,
      bamLibraryDesigns=bamLibraryDesigns,
      svFiles=svFiles,
      workflowNames=workflowNames,
      svLibraryDesigns=svLibraryDesigns
  }

  ## Feed output of generateConfig to input of runMavis
  File mavisConfig = generateConfig.jsonFile

  call runMavis { 
    input: 
      configFile = mavisConfig,
      outputFileNamePrefix = sanitizedSID,
      modules = mavis_modules
  }


  meta {
    author: "Hannah Driver"
    email: "hdriver@oicr.on.ca"
    description: "MAVIS workflow, annotation of structural variants. An application framework for the rapid generation of structural variant consensus, able to visualize the genetic impact and context as well as process both genome and transcriptome data."
    dependencies: [
      {
        name: "mavis/3.1.0",
        url: "http://mavis.bcgsc.ca/"
      }
    ]
    output_meta: {
      summary: "File with copy number variants, native varscan format",
      drawings: "Plots generated with MAVIS, collected into a single tar.gz archive",
      nscvWT: "Whole transcriptome non-synonymous coding variants. The output file is only generated if variants are found",
      nscvWG: "Whole genome non-synonymous coding variants. The output file is only generated if variants are found"
    }
  }


  output {
    File summary = runMavis.summary
    File drawings = runMavis.drawings
    File? nscvWT   = runMavis.nscvWT
    File? nscvWG   = runMavis.nscvWG
  }

}

task filterDelly {
  input {
    String outputFileNamePrefix
    Array[String] svFiles
    Array[String] workflowNames
    Array[String] svLibraryDesigns
    String modules
    String megafusionExecutable
    String megafusionArriba
    String megafusionStarfusion
    String survivorExecutable
    Int maxLines = 2000
    Int variantSupport = 10
    Int jobMemory = 24
    Int timeout = 6
  }

  parameter_meta {
    outputFileNamePrefix: "Sample ID, this is provided to mavis and cannot include reserved characters [;,_\\s] "
    svFiles: "Array of SV calls"
    workflowNames: "List of workflow names for SV inputs (e.g. delly, starfusion, arriba)"
    svLibraryDesigns: "List of library designs (e.g. WG, WT)"
    modules: "Modules needed to filter delly file"
    megafusionExecutable: "Path to MegaFusion executable"
    megafusionArriba: "Path to MegaFusion arriba.json file"
    megafusionStarfusion: "Path to MegaFusion starfusion.json file"
    survivorExecutable: "Path to SURVIVOR executable"
    maxLines: "Maximum number of lines a delly file can have before needing filtering. Default is 2000"
    variantSupport: "Paired-end support for structural variants, in pairs. Default is 10"
    jobMemory: "Memory allocated for this job"
    timeout: "Timeout in hours, needed to override imposed limits"
  }

  command <<<

    #See GRD-744 for breakdown of this task

    set -eu -o pipefail

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

            if lines > ~{maxLines}:
                #Check if other SV callers exist or else survivor can't be run
                if len(svFiles) > 1:
                    #Run megafusion
                    for index, name in enumerate(workflowNames):
                        if name.lower() == "arriba":
                            arriba_command = f'python3 ~{megafusionExecutable} --json ~{megafusionArriba} --fusion {svFiles_escaped[index]} > arriba.vcf'
                            subprocess.run(arriba_command, shell=True)
                        if name.lower() == "starfusion":
                            starfusion_command = f'python3 ~{megafusionExecutable} --json ~{megafusionStarfusion} --fusion {svFiles_escaped[index]} > starfusion.vcf'
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
                    survivor_command = f'"~{survivorExecutable}" merge "survivor_input.txt" 1000 2 0 0 0 1 merged.vcf'
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
                subprocess.run(['bcftools', 'view', '-i', f'FILTER="PASS" & INFO/PE>~{variantSupport}', '-O', 'z', original_delly, '-o', 'filtered_delly.vcf.gz'])

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
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File mavisDelly = '~{outputFileNamePrefix}_mavis_delly.vcf'
  }
}


task generateConfig {
  input {
    String outputFileNamePrefix
    String diseaseStatus
    String alignerReference
    String annotations
    String dgvAnnotation
    String masking
    String referenceGenome
    String templateMetadata
    String reference
    String modules
    Array[String] bams
    Array[String] bamLibraryDesigns
    Array[String] svFiles
    Array[String] workflowNames
    Array[String] svLibraryDesigns
    Boolean? drawFusionsOnly
    Int? minClustersPerFile
    Boolean? uninformativeFilter
    Int? filterMinFlankingReads
    Int? filterMinLinkingSplitReads
    Int? filterMinRemappedReads
    Int? filterMinSpanningReads
    Boolean? filterTransHomopolymers
    Int jobMemory = 6
    Int timeout = 6
  }


  parameter_meta {
    outputFileNamePrefix: "Sample ID, this is provided to maivs and cannot include reseerved characters [;,_\\s] "
    diseaseStatus: "Tissue status. For example: diseased"
    alignerReference: "References in 2bit (compressed) format, used by MAVIS aligner"
    annotations: ".json file with annotations for MAVIS"
    dgvAnnotation: "The DGV annotations help to deal with variants found in normal tissue"
    masking: "Masking data in .tab format"
    referenceGenome: "Path to fasta file with genomic assembly"
    templateMetadata: "Chromosome Band Information, used for visualization"
    reference: "The genome reference build. for example: hg19, hg38"
    modules: "Modules needed to run MAVIS"
    bams: "Array of input bam file paths"
    bamLibraryDesigns: "List of library designs (e.g. WG, WT)"
    svFiles: "Array of SV calls"
    workflowNames: "List of workflow names for SV iputs (e.g. delly, starfusion, arriba)"
    svLibraryDesigns: "List of library designs (e.g. WG, WT)"
    drawFusionsOnly: "flag for MAVIS visualization control" 
    minClustersPerFile: "Determines the way parallel calculations are organized"
    uninformativeFilter: "If enabled, only interested in events inside genes, speeds up calculations"
    filterMinFlankingReads: "Minimum number of flanking pairs for a call by flanking pairs"
    filterMinLinkingSplitReads: "Minimum number of linking split reads for a call by split reads"
    filterMinRemappedReads: "Minimum number of remapped reads for a call by contig"
    filterMinSpanningReads: "Minimum number of spanning reads for a call by spanning reads"
    filterTransHomopolymers: "When enabled, transcript sequences containing homopolymer regions are removed"
    jobMemory: "Memory allocated for this job"
    timeout: "Timeout in hours, needed to override imposed limits"
  }


  ## Declare default values for config file flags
  Boolean drawFusionsActual = select_first([drawFusionsOnly, false])
  Int minClustersActual = select_first([minClustersPerFile, 100])
  Boolean uninformativeFilterActual = select_first([uninformativeFilter, true])
  Int filterMinFlankingActual = select_first([filterMinFlankingReads, 10])
  Int filterMinLinkingActual = select_first([filterMinLinkingSplitReads, 1])
  Int filterMinRemappedActual = select_first([filterMinRemappedReads, 5])
  Int filterMinSpanningActual = select_first([filterMinSpanningReads, 5])
  Boolean filterTransActual = select_first([filterTransHomopolymers, false])



  command <<<

    set -eu -o pipefail

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
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File jsonFile = "config.json"
  }
}


task runMavis {
  input {
    File configFile
    String outputFileNamePrefix
    String modules
    Int jobMemory = 96
    Int timeout = 24
  }


  parameter_meta {
    configFile: ".json configuration file for mavis"
    outputFileNamePrefix: "sample ID, this is provided to maivs and cannot include reseerved characters [;,_\\s] "
    modules: "modules needed to run MAVIS"
    jobMemory: "Memory allocated for this job"
    timeout: "Timeout in hours, needed to override imposed limits"
  }


  command <<<

    set -eu -o pipefail
    
    snakemake --jobs 40 --configfile=~{configFile} -s $MAVIS_ROOT/bin/Snakefile


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

  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File drawings  = "~{outputFileNamePrefix}.mavis_drawings.zip"
    File summary   = "~{outputFileNamePrefix}.mavis_summary.tab"
    File? nscvWT   = "~{outputFileNamePrefix}.WT_non-synonymous_coding_variants.tab"
    File? nscvWG   = "~{outputFileNamePrefix}.WG_non-synonymous_coding_variants.tab"
  }

}
