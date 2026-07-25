# mavis3

MAVIS workflow, annotation of structural variants. An application framework for the rapid generation of structural variant consensus, able to visualize the genetic impact and context as well as process both genome and transcriptome data.

## Overview

## Dependencies

* [mavis 3.1.2](http://mavis.bcgsc.ca/)


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
`filterDelly.maxLines`|Int|2000|Maximum number of lines a delly file can have before needing filtering. Default is 2000
`filterDelly.variantSupport`|Int|10|Paired-end support for structural variants, in pairs. Default is 10
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
`setup.containerCommand`|String|"singularity exec -B /.mounts/ -B /scratch3/ $MAVIS_ROOT/bin/mavis.sif"|Command used to enter the mavis container. The mavis module ships only a snakemake wrapper, so mavis and its bundled python helpers are invoked through singularity directly
`setup.singularityTmpDir`|String|"/tmp"|Directory singularity extracts the image into. Must be node-local, an NFS-backed default costs minutes per task
`setup.maxBins`|Int|100000|Maximum value for the sample_bin_size parameter if the config fails to build, Default is 100000
`setup.jobMemory`|Int|16|Memory allocated for this job
`setup.cores`|Int|1|Number of cores allocated for this job
`setup.timeout`|Int|12|Timeout in hours, needed to override imposed limits
`convert.containerCommand`|String|"singularity exec -B /.mounts/ -B /scratch3/ $MAVIS_ROOT/bin/mavis.sif"|Command used to enter the mavis container. The mavis module ships only a snakemake wrapper, so mavis and its bundled python helpers are invoked through singularity directly
`convert.singularityTmpDir`|String|"/tmp"|Directory singularity extracts the image into. Must be node-local, an NFS-backed default costs minutes per task
`convert.jobMemory`|Int|16|Memory allocated for this job
`convert.cores`|Int|1|Number of cores allocated for this job
`convert.timeout`|Int|6|Timeout in hours, needed to override imposed limits
`cluster.containerCommand`|String|"singularity exec -B /.mounts/ -B /scratch3/ $MAVIS_ROOT/bin/mavis.sif"|Command used to enter the mavis container. The mavis module ships only a snakemake wrapper, so mavis and its bundled python helpers are invoked through singularity directly
`cluster.singularityTmpDir`|String|"/tmp"|Directory singularity extracts the image into. Must be node-local, an NFS-backed default costs minutes per task
`cluster.jobMemory`|Int|16|Memory allocated for this job
`cluster.cores`|Int|1|Number of cores allocated for this job
`cluster.timeout`|Int|6|Timeout in hours, needed to override imposed limits
`validateAndAnnotate.batchName`|String|basename(batchFile,".tab")|Name of the batch, taken from the batch file name and used for the output directories
`validateAndAnnotate.containerCommand`|String|"singularity exec -B /.mounts/ -B /scratch3/ $MAVIS_ROOT/bin/mavis.sif"|Command used to enter the mavis container. The mavis module ships only a snakemake wrapper, so mavis and its bundled python helpers are invoked through singularity directly
`validateAndAnnotate.singularityTmpDir`|String|"/tmp"|Directory singularity extracts the image into. Must be node-local, an NFS-backed default costs minutes per task
`validateAndAnnotate.jobMemory`|Int|18|Memory allocated for this job
`validateAndAnnotate.cores`|Int|2|Number of cores allocated for this job
`validateAndAnnotate.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`pairing.containerCommand`|String|"singularity exec -B /.mounts/ -B /scratch3/ $MAVIS_ROOT/bin/mavis.sif"|Command used to enter the mavis container. The mavis module ships only a snakemake wrapper, so mavis and its bundled python helpers are invoked through singularity directly
`pairing.singularityTmpDir`|String|"/tmp"|Directory singularity extracts the image into. Must be node-local, an NFS-backed default costs minutes per task
`pairing.jobMemory`|Int|16|Memory allocated for this job
`pairing.cores`|Int|1|Number of cores allocated for this job
`pairing.timeout`|Int|12|Timeout in hours, needed to override imposed limits
`mavisSummary.containerCommand`|String|"singularity exec -B /.mounts/ -B /scratch3/ $MAVIS_ROOT/bin/mavis.sif"|Command used to enter the mavis container. The mavis module ships only a snakemake wrapper, so mavis and its bundled python helpers are invoked through singularity directly
`mavisSummary.singularityTmpDir`|String|"/tmp"|Directory singularity extracts the image into. Must be node-local, an NFS-backed default costs minutes per task
`mavisSummary.jobMemory`|Int|16|Memory allocated for this job
`mavisSummary.cores`|Int|1|Number of cores allocated for this job
`mavisSummary.timeout`|Int|12|Timeout in hours, needed to override imposed limits


### Outputs

Output | Type | Description | Labels
---|---|---|---
`summary`|File|File with copy number variants, native varscan format|vidarr_label: summary
`drawings`|File|Plots generated with MAVIS, collected into a single tar.gz archive|vidarr_label: drawings
`nscvWT`|File?|Whole transcriptome non-synonymous coding variants. The output file is only generated if variants are found|vidarr_label: nscvWT
`nscvWG`|File?|Whole genome non-synonymous coding variants. The output file is only generated if variants are found|vidarr_label: nscvWG


## Commands
This section lists command(s) run by mavis3 workflow

* Running mavis3

```

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
```
```

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
        mavisRoot = str(os.environ['HG38V110_MAVIS_ROOT'])

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
            "bam_stats.sample_bin_size": 1000,
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

            if name.lower() == "gridss":
                if "gridss" not in jsonDict["convert"]:
                    jsonDict["convert"]["gridss"] = {
                        "assume_no_untemplated": True,
                        "file_type": "vcf",
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
            if name.lower() == "gridss":
                jsonDict["libraries"]["WG." + "~{outputFileNamePrefix}"]["assign"].append("gridss")
            if name.lower() == "starfusion":
                jsonDict["libraries"]["WT." + "~{outputFileNamePrefix}"]["assign"].append("starfusion")
            if name.lower() == "arriba":
                jsonDict["libraries"]["WT." + "~{outputFileNamePrefix}"]["assign"].append("arriba")

        with open("config.json", 'w') as jsonFile:
            json.dump(jsonDict, jsonFile)

    CODE
```
```

    set -euo pipefail

    ## Singularity on this cluster cannot loop-mount the SIF, so every invocation extracts
    ## the image. Keeping that on node-local disk rather than the NFS-backed default TMPDIR
    ## takes it from ~140s to ~3s per task.
    export SINGULARITY_TMPDIR="~{singularityTmpDir}"

    cp ~{configFile} config.raw.json

    ## mavis setup samples the input bams to add bam_stats and fills in schema defaults.
    ## If the requested bin size cannot be satisfied, retry once with a larger one.
    if ! ~{containerCommand} mavis setup --config config.raw.json --outputfile config.initialized.json; then
      sed -i 's/bin_size": 1000/bin_size": ~{maxBins}/' config.raw.json
      ~{containerCommand} mavis setup --config config.raw.json --outputfile config.initialized.json
    fi

    ## Runs inside the container because it imports the mavis_config helpers shipped with
    ## mavis itself, rather than reimplementing them.
    ~{containerCommand} python3 <<CODE
    import json

    from mavis_config import get_library_inputs, guess_total_batches

    with open("config.initialized.json") as handle:
        config = json.load(handle)

    #Every stage below runs in its own Cromwell directory and is given an explicit --output,
    #so output_dir must stay relative rather than pointing at this task's directory.
    config["output_dir"] = "output_dir_full"

    #mavis cluster requires total_batches per library: it is what decides how many batch
    #files the library is split into. The upstream Snakefile computes this at graph-build
    #time, so it is neither written by generateConfig nor added by mavis setup.
    for library in config["libraries"]:
        if "total_batches" in config["libraries"][library]:
            continue
        config["libraries"][library]["total_batches"] = guess_total_batches(
            config, get_library_inputs(config, library)
        )

    with open("config.initialized.json", "w") as handle:
        json.dump(config, handle, sort_keys=True, indent="  ")

    #The library and converter alias lists drive the scatters in the workflow body. Taking
    #them from the initialized config keeps them consistent with what mavis will accept.
    with open("libraries.txt", "w") as handle:
        for library in sorted(config["libraries"]):
            handle.write(library + "\n")

    with open("aliases.txt", "w") as handle:
        for alias in sorted(config.get("convert", {})):
            handle.write(alias + "\n")
    CODE
```
```

    set -euo pipefail

    export SINGULARITY_TMPDIR="~{singularityTmpDir}"

    mkdir -p converted_outputs output_dir_full

    ## The conversion flags are read from the initialized config rather than passed in:
    ## mavis setup is what fills in the schema defaults, and strand_specific in particular
    ## is never written by generateConfig.
    python3 <<CODE
    import json
    import shlex

    with open("~{configFile}") as handle:
        settings = json.load(handle)["convert"]["~{converterAlias}"]

    arguments = [
        "~{containerCommand} mavis",
        "convert",
        "--file_type", str(settings["file_type"]),
        "--strand_specific", str(settings["strand_specific"]),
        "--assume_no_untemplated", str(settings["assume_no_untemplated"]),
        "--inputs",
    ]
    arguments.extend(shlex.quote(str(path)) for path in settings["inputs"])
    arguments.extend(["--outputfile", "converted_outputs/~{converterAlias}.tab"])

    with open("convert.sh", "w") as handle:
        handle.write("#!/bin/bash\n")
        handle.write("set -euo pipefail\n")
        handle.write(" ".join(arguments) + "\n")
    CODE

    bash convert.sh
```
```

    set -euo pipefail

    export SINGULARITY_TMPDIR="~{singularityTmpDir}"

    mkdir -p output_dir_full

    ## Each entry of the library's assign list names a converter whose output is called
    ## <alias>.tab; anything not produced by a converter is already a path. This mirrors
    ## get_cluster_inputs in the MAVIS Snakefile.
    python3 <<CODE
    import json
    import os
    import shlex

    with open("~{configFile}") as handle:
        config = json.load(handle)

    converted = {}
    for path in "~{sep=' ' convertedFiles}".split():
        converted[os.path.basename(path)] = path

    inputs = []
    for alias in config["libraries"]["~{library}"]["assign"]:
        inputs.append(converted.get(alias + ".tab", alias))

    arguments = [
        "~{containerCommand} mavis",
        "cluster",
        "--config", shlex.quote("~{configFile}"),
        "--library", shlex.quote("~{library}"),
        "--inputs",
    ]
    arguments.extend(shlex.quote(path) for path in inputs)
    arguments.extend(["--output", shlex.quote("cluster")])

    with open("cluster.sh", "w") as handle:
        handle.write("#!/bin/bash\n")
        handle.write("set -euo pipefail\n")
        handle.write(" ".join(arguments) + "\n")
    CODE

    bash cluster.sh
```
```

    set -euo pipefail

    export SINGULARITY_TMPDIR="~{singularityTmpDir}"

    mkdir -p output_dir_full

    ## validate and annotate are one-to-one per batch, so they are fused into a single task:
    ## splitting them would localize the same intermediates twice for no added tracking.
    ~{containerCommand} mavis validate \
      --config ~{configFile} \
      --library ~{library} \
      --inputs ~{batchFile} \
      --output validate/~{batchName}

    ~{containerCommand} mavis annotate \
      --config ~{configFile} \
      --library ~{library} \
      --inputs validate/~{batchName}/validation-passed.tab \
      --output annotate/~{batchName}

    if [ ! -f annotate/~{batchName}/MAVIS.COMPLETE ]; then
      echo "ERROR: annotate did not complete for ~{library} ~{batchName}" >&2
      exit 1
    fi
```
```

    set -euo pipefail

    export SINGULARITY_TMPDIR="~{singularityTmpDir}"

    ## mavis pairing writes straight into --output without creating it: snakemake used to
    ## create the parent directories of a rule's declared outputs beforehand.
    mkdir -p output_dir_full pairing

    ~{containerCommand} mavis pairing \
      --config ~{configFile} \
      --inputs ~{sep=' ' annotations} \
      --output pairing

    if [ ! -f pairing/MAVIS.COMPLETE ]; then
      echo "ERROR: pairing did not complete" >&2
      exit 1
    fi
```
```

    set -euo pipefail

    export SINGULARITY_TMPDIR="~{singularityTmpDir}"

    ## As with pairing, create the output directory rather than relying on the caller to
    ## have made it.
    mkdir -p output_dir_full summary

    ~{containerCommand} mavis summary \
      --config ~{configFile} \
      --inputs ~{pairedFile} \
      --output summary

    if [ ! -f summary/MAVIS.COMPLETE ]; then
      echo "MAVIS job finished but THERE ARE NO RESULTS" >&2
      exit 1
    fi

    ### create an empty zip file, which will be updated with drawings and legends.  if there are none, than the empty file is provisioned out
    echo | zip -q > ~{outputFileNamePrefix}.mavis_drawings.zip && zip -dq ~{outputFileNamePrefix}.mavis_drawings.zip -

    ### the leading "" keeps the loop valid when no drawings were produced at all
    : > drawings.list
    for drawing in "" ~{sep=' ' drawingFiles}; do
      if [ -n "$drawing" ]; then
        printf '%s\n' "$drawing" >> drawings.list
      fi
    done

    ### add every drawing to the archive, flattened, as the single-task version did
    if [ -s drawings.list ]; then
      zip -qj ~{outputFileNamePrefix}.mavis_drawings.zip -@ < drawings.list
    fi

    ### there should be a single mavis_summary_all files
    cp summary/mavis_summary_all_*.tab ~{outputFileNamePrefix}.mavis_summary.tab

    ### non-synonymous coding variants are separate into WG or WT files; each may or may not be produced
    if [ -e summary/mavis_summary_WG.*_non-synonymous_coding_variants.tab ];then
      cp summary/mavis_summary_WG.*_non-synonymous_coding_variants.tab ~{outputFileNamePrefix}.WG_non-synonymous_coding_variants.tab
    fi
    if [ -e summary/mavis_summary_WT.*_non-synonymous_coding_variants.tab ];then
      cp summary/mavis_summary_WT.*_non-synonymous_coding_variants.tab ~{outputFileNamePrefix}.WT_non-synonymous_coding_variants.tab
    fi
```
## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
