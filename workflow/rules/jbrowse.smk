# ─────────────────────────────────────────────────────────────────────────────
# 1. Add Assembly
# ─────────────────────────────────────────────────────────────────────────────

rule jbrowse_add_assembly:
    output:
        config=temp("results/jbrowse/config_{genome}_assembly.json"),
    log:
        "results/jbrowse/add_assembly_{genome}.log",
    conda:
        "../envs/jbrowse.yml"
    resources:
        file_lock=1,
    params:
        url_prefix=lambda wc: config["genome_dir"],
        extra=lambda wc: config["jbrowse"]["add_assembly"][wc.genome],
    message:
        "add {wildcards.genome} assemblies to jbrowse"
    shell:
        """
        (
            # Add to jbrowse
            url="{params.url_prefix}{wildcards.genome}.2bit"
            jbrowse add-assembly $url \
                --type twoBit \
                --target {output.config} \
                ${params.extra}
        ) >{log} 2>&1
        """


# ─────────────────────────────────────────────────────────────────────────────
# 2. Add Annotation
# ─────────────────────────────────────────────────────────────────────────────

rule jbrowse_add_anno:
    input:
        config="results/jbrowse/config_{genome}_assembly.json",
    output:
        config=temp("results/jbrowse/config_{genome}_anno.json"),
    log:
        "results/jbrowse/add_anno_{genome}.log",
    conda:
        "../envs/jbrowse.yml"
    resources:
        file_lock=1,
    params:
        url_prefix=lambda wc: config["genome_dir"],
        extra=lambda wc: config["jbrowse"]["add_anno"][wc.genome],
    message:
        "add {wildcards.genome} annotations to jbrowse"
    shell:
        """
        (

            # Get url/path to annotation files
            path_to_gff="{params.url_prefix}{wildcards.genome}.sorted.gff.gz"
            path_to_tbi="{params.url_prefix}{wildcards.genome}.sorted.gff.gz.tbi"

            # add to jbrowse
            cp {input.config} {output.config}
            jbrowse add-track "$path_to_gff" \
                --indexFile "$path_to_tbi" \
                --target {output.config} \
                --assemblyNames {wildcards.genome} \
                {params.extra}
        ) >{log} 2>&1
        """


# ─────────────────────────────────────────────────────────────────────────────
# 3. Add BigWigs
# ─────────────────────────────────────────────────────────────────────────────


rule jbrowse_add_bw:
    input:
        config="results/jbrowse/config_{genome}_anno.json",
    output:
        config=temp("results/jbrowse/config_{genome}_{study}_bw.json"),
    log:
        "results/jbrowse/add_bw_{genome}_{study}.log",
    conda:
        "../envs/jbrowse.yml"
    resources:
        file_lock=1,
    params:
        plus_bw=lambda wc: expand(
            "results/deeptools/coverage/{sample}.{genome}.plus.bw",
            sample=samples[samples.study == wc.study]["secondary_name"],
            genome=wc.genome
        ),
        minus_bw=lambda wc: expand(
            "results/deeptools/coverage/{sample}.{genome}.minus.bw",
            sample=samples[samples.study == wc.study]["secondary_name"],
            genome=wc.genome
        ),
        display_names=lambda wc: list(samples[samples.study == wc.study]["sample"]),
        url_prefix=lambda wc: config["sample_urls"][wc.study],
        color=lambda wc: config["color"][wc.study],
        extra=lambda wc: config["jbrowse"]["add_bw"][wc.study],
    message:
        "add BigWigs from {wildcards.study} to {wildcards.genome} config.json"
    shell:
        """
        (
            cp {input.config} {output.config}

            display_names=({params.display_names})
            plus_bws=({params.plus_bw})
            minus_bws=({params.minus_bw})

            length=${{#display_names[@]}}

            for ((i = 0; i < length; i++)); do

                samplename="${{display_names[$i]}}"
                plus_bw="${{plus_bws[$i]}}"
                minus_bw="${{minus_bws[$i]}}"

                path_to_plus_bw="{params.url_prefix}/${{plus_bw}}"
                path_to_minus_bw="{params.url_prefix}/${{minus_bw}}"

                jbrowse add-track $path_to_plus_bw \
                    --target {output.config} \
                    --name ${{samplename}}_plus_bw \
                    --trackId ${{samplename}}_{wildcards.genome}_plus_bw \
                    --assemblyNames {wildcards.genome} \
                    --config '{{"displays":[{{"type":"LinearWiggleDisplay","color":"{params.color}"}}]}}' \
                    {params.extra}

                jbrowse add-track $path_to_minus_bw \
                    --target {output.config} \
                    --name ${{samplename}}_minus_bw \
                    --trackId ${{samplename}}_{wildcards.genome}_minus_bw \
                    --assemblyNames {wildcards.genome} \
                    --config '{{"displays":[{{"type":"LinearWiggleDisplay","inverted":true,"color":"{params.color}"}}]}}' \
                    {params.extra}
            done

        ) >{log} 2>&1
        """


# ─────────────────────────────────────────────────────────────────────────────
# 4. Add cram files
# ─────────────────────────────────────────────────────────────────────────────


rule jbrowse_add_cram:
    input:
        config="results/jbrowse/config_{genome}_{study}_bw.json",
    output:
        config=temp("results/jbrowse/config_{genome}_{study}_cram.json"),
    log:
        "results/jbrowse/add_cram_{genome}_{study}.log",
    conda:
        "../envs/jbrowse.yml"
    resources:
        file_lock=1,
    params:
        cram=lambda wc: expand(
            "results/processed_alignment/cram/{sample}_{genome}.cram",
            sample=samples[samples.study == wc.study]["secondary_name"],
            genome=wc.genome
        ),
        display_names=lambda wc: list(samples[samples.study == wc.study]["sample"]),
        url_prefix=lambda wc: config["sample_urls"][wc.study],
        extra=lambda wc: config["jbrowse"]["add_cram"][wc.study],
    message:
        "add plus cram tracks to {wildcards.genome} config.json"
    shell:
        """
        (            
            cp {input.config} {output.config}
            display_names=({params.display_names})
            cram_files=({params.cram})
            length=${{#display_names[@]}}

            for ((i = 0; i < length; i++)); do

                samplename="${{display_names[$i]}}"
                cram="${{cram_files[$i]}}"
                path_to_cram="{params.url_prefix}/${{cram}}"

                # Add to jbrowse
                jbrowse add-track $path_to_cram \
                    --indexFile $i.crai \
                    --target {output.config} \
                    --name ${{samplename}}_cram \
                    --trackId ${{samplename}}_{wildcards.genome}_cram \
                    --assemblyNames {wildcards.genome} \
                    --config '{{"displays":[{{"type":"LinearPileupDisplay", "showLegend": true, "colorBySetting": {{"type": "stranded"}}}}]}}' \
                    {params.extra}
            done
        ) >{log} 2>&1
        """


# ─────────────────────────────────────────────────────────────────────────────
# 5. Merge config.jsons
# ─────────────────────────────────────────────────────────────────────────────


rule jbrowse_merge_jsons:
    input:
        config_files=expand(
            "results/jbrowse/config_{genome}_{study}_cram.json",
            genome=list(config["jbrowse"]["add_assembly"].keys()),
            study=list(config["sample_urls"].keys()),
        ),
    output:
        "results/jbrowse/config.json",
    log:
        "results/jbrowse/merge_config_json.log",
    conda:
        "../envs/jbrowse.yml"
    message:
        "merge config files"
    shell:
        """
        jq -s '
            {{
                assemblies:                  map(.assemblies)                  | add | unique_by(.name),
                tracks:                      map(.tracks)                      | add | unique_by(.trackId),
                connections:                 map(.connections // [])           | add | unique_by(.name),
                aggregateTextSearchAdapters: map(.aggregateTextSearchAdapters // []) | add | unique_by(.textSearchAdapterId),
                configuration:               map(.configuration)               | add,
                defaultSession:              .[0].defaultSession
            }}
        ' {input.config_files} >{output}
        """