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

  ## Setting modules and reference paths by genome
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
    String svFiles = "~{s.svFile}"
    String workflowNames = s.workflowName
    String svLibraryDesigns = s.libraryDesign
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
    email: "hannah.driver@oicr.on.ca"
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
            "convert": {}
        }

        for index, name in enumerate(workflowNames):
            if name.lower() == "delly":
                entry = {
                    "delly": {
                        "assume_no_untemplated": True,
                        "file_type": "delly",
                        "inputs": [
                            str(svFiles[index])
                        ]
                    }
                }
                jsonDict["convert"].append(entry)

            if name.lower() == "starfusion":
                entry = {
                    "starfusion": {
                        "assume_no_untemplated": True,
                        "file_type": "starfusion",
                        "inputs": [
                            str(svFiles[index])
                        ]
                    }
                }
                jsonDict["convert"].append(entry)

            if name.lower() == "arriba":
                entry = {
                    "arriba": {
                        "assume_no_untemplated": True,
                        "file_type": "arriba",
                        "inputs": [
                            str(svFiles[index])
                        ]
                    }
                }
                jsonDict["convert"].append(entry)
  
        for index, bam in enumerate(bams):
            if bamLibraryDesigns[index] == "WG":
                jsonDict["libraries"] = {
                    "WG." + "~{outputFileNamePrefix}": {
                        "assign": [
                            "delly"
                        ],
                        "bam_file": bams[index],
                        "disease_status": "~{diseaseStatus}",
                        "protocol": "genome"
                    }
                }
            if bamLibraryDesigns[index] == "WT":
                jsonDict["libraries"] = {
                    "WT." + "~{outputFileNamePrefix}": {
                        "assign": [],
                        "bam_file": bams[index],
                        "disease_status": "~{diseaseStatus}",
                        "protocol": "transcriptome",
                        "strand_specific": True
                    }
                }

        for name in workflowNames:
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
    Int jobMemory = 36 
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
    
    snakemake --jobs 100 --configfile=~{configFile} -s $MAVIS_ROOT/bin/Snakefile


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
