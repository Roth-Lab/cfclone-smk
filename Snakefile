import os
from snakemake.io import expand, directory, temp
from snakemake.utils import min_version

min_version("9.16")


conda: "envs/global.yaml"


include: "utils.py"


config = ConfigManager(config)


rule all:
    localrule: True
    input:
        config.pipeline_files,


rule run_cfclone:
    input:
        c=config.clone_cn_file,
        i=config.ctdna_file,
    output:
        e=directory(config.exec_dir_template),
        f=config.fit_template,
    params:
        c=config.num_chains,
        r=config.num_rounds,
        rt=config.get_cfclone_run_type_args,
    benchmark:
        config.get_benchmark_file(config.fit_template)
    conda:
        "envs/cfclone.yaml"
    log:
        config.get_log_file(config.fit_template),
    threads: config.num_threads
    shell:
        "(cfclone fit "
        "-c {input.c} "
        "-i {input.i} "
        "-o {output.f} "
        "-t {threads} "
        "--exec-dir {output.e} "
        "--num-chains {params.c} "
        "--num-rounds {params.r} "
        "--seed {wildcards.seed} "
        "{params.rt})  >{log} 2>&1"


rule write_ancestral_prevalances:
    input:
        i=lambda wildcards: expand(config.fit_template, run_type="full", seed=wildcards.seed),
        t=config.clone_tree_file,
    output:
        o=config.ancestral_prevalence_file,
        t=config.clone_prevalence_tree_file,
    conda:
        "envs/cfclone.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.ancestral_prevalence_file),
    shell:
        "(cfclone write-ancestral-prevalences -c {input.t} -i {input.i} -o {output.o} -t {output.t}) >{log} 2>&1"


rule write_dominance_prob:
    input:
        fit=lambda wildcards: expand(config.fit_template, run_type="full", seed=wildcards.seed),
    output:
        config.dominance_prob_file,
    conda:
        "envs/cfclone.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.dominance_prob_file),
    shell:
        "(cfclone write-dominance-prob -i {input.fit} -o {output}) >{log} 2>&1"


rule write_pairwise_ranks_file:
    input:
        fit=lambda wildcards: expand(config.fit_template, run_type="full", seed=wildcards.seed),
    output:
        config.pairwise_ranks_file,
    conda:
        "envs/cfclone.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.pairwise_ranks_file),
    shell:
        "(cfclone write-pairwise-ranks -i {input.fit} -o {output}) >{log} 2>&1"


rule write_parameter_summaries_file:
    input:
        fit=lambda wildcards: expand(config.fit_template, run_type="full", seed=wildcards.seed),
    output:
        config.parameter_summaries_file,
    conda:
        "envs/cfclone.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.parameter_summaries_file),
    shell:
        "(cfclone write-parameter-summaries -i {input.fit} -o {output}) >{log} 2>&1"


rule write_summary_file:
    input:
        fit=lambda wildcards: expand(config.fit_template, run_type="full", seed=wildcards.seed),
    output:
        config.summary_file,
    conda:
        "envs/cfclone.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.summary_file),
    shell:
        "(cfclone write-summary -i {input.fit} -o {output}) >{log} 2>&1"


rule write_tumour_content_file:
    input:
        fit=lambda wildcards: expand(config.fit_template, run_type="full", seed=wildcards.seed),
    output:
        config.tumour_content_file,
    conda:
        "envs/cfclone.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.tumour_content_file),
    shell:
        "(cfclone write-tumour-content -i {input.fit} -o {output}) >{log} 2>&1"


rule write_evidence:
    input:
        config.fit_template,
    output:
        temp(config.run_type_evidence_template),
    conda:
        "envs/cfclone.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.run_type_evidence_template),
    shell:
        "(echo -n {wildcards.run_type}, > {output}; "
        "cfclone print-model-evidence -i{input} >> {output}) 2>{log}"


rule merge_evidence:
    input:
        i=lambda wildcards: expand(config.run_type_evidence_template, run_type=["full", "normal"], seed=wildcards.seed),
        script=workflow.source_path("scripts/merge_evidence.py"),
    output:
        config.evidence_file,
    conda:
        "envs/python.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.evidence_file),
    shell:
        "(python {input.script} -i {input.i} -o {output}) >{log} 2>&1"


rule plot_clone_prevalences:
    input:
        d=config.ancestral_prevalence_file,
        t=config.clone_prevalence_tree_file,
        script=workflow.source_path("scripts/plot_clone_prevalences.py"),
    output:
        config.clone_prevalences_plot,
    conda:
        "envs/plot.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.clone_prevalence_tree_file),
    shell:
        "(python {input.script} -d {input.d} -t {input.t} -o {output}) >{log} 2>&1"


rule plot_fit:
    input:
        i=config.parameter_summaries_file,
        script=workflow.source_path("scripts/plot_fit.py"),
    output:
        config.fit_plot,
    conda:
        "envs/plot.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.fit_plot),
    shell:
        "(python {input.script} -i {input.i} -o {output}) >{log} 2>&1"


rule plot_pairwsie_ranks:
    input:
        i=config.pairwise_ranks_file,
        t=config.clone_tree_file,
        script=workflow.source_path("scripts/plot_pairwise_ranks.py"),
    output:
        config.pairwise_ranks_plot,
    conda:
        "envs/plot.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.pairwise_ranks_plot),
    shell:
        "(python {input.script} -i {input.i} -t {input.t} -o {output}) >{log} 2>&1"


rule merge_evidence_restarts:
    input:
        i=expand(config.evidence_file, seed=config.seeds),
        script=workflow.source_path("scripts/merge_restart_tables.py"),
    output:
        config.merged_evidence_file,
    params:
        seeds=" ".join([str(x) for x in config.seeds]),
    conda:
        "envs/python.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.merged_evidence_file),
    shell:
        "(python {input.script} -i {input.i} -o {output} -s {params.seeds}) >{log} 2>&1"


rule merge_prevalence_restarts:
    input:
        i=expand(config.ancestral_prevalence_file, seed=config.seeds),
        script=workflow.source_path("scripts/merge_restart_tables.py"),
    output:
        config.merged_prevalence_file,
    params:
        seeds=" ".join([str(x) for x in config.seeds]),
    conda:
        "envs/python.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.merged_prevalence_file),
    shell:
        "(python {input.script} -i {input.i} -o {output} -s {params.seeds}) >{log} 2>&1"


rule merge_summary_restarts:
    input:
        i=expand(config.summary_file, seed=config.seeds),
        script=workflow.source_path("scripts/merge_restart_tables.py"),
    output:
        config.merged_summary_file,
    params:
        seeds=" ".join([str(x) for x in config.seeds]),
    conda:
        "envs/python.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.merged_summary_file),
    shell:
        "(python {input.script} -i {input.i} -o {output} -s {params.seeds}) >{log} 2>&1"


rule merge_tumour_content_restarts:
    input:
        i=expand(config.tumour_content_file, seed=config.seeds),
        script=workflow.source_path("scripts/merge_restart_tables.py"),
    output:
        config.merged_tumour_content_file,
    params:
        seeds=" ".join([str(x) for x in config.seeds]),
    conda:
        "envs/python.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.merged_tumour_content_file),
    shell:
        "(python {input.script} -i {input.i} -o {output} -s {params.seeds}) >{log} 2>&1"


rule save_run_configuration:
    localrule: True
    output:
        config.experiment_configuration,
    params:
        config.config_yaml,
    log:
        config.log_dir.joinpath("save_run_configuration.log"),
    conda:
        "envs/python.yaml"
    shell:
        "(echo \"{params}\" > {output}) 2>{log}"
