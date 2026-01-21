import os

from snakemake.utils import min_version, validate

min_version("9.4")


conda: workflow.source_path("envs/global.yaml")


validate(config, "schemas/config.schema.yaml")

# from utils import ConfigManager


include: workflow.source_path("utils.py")


config = ConfigManager(config)


rule all:
    input:
        config.pipeline_files,
    localrule: True


if config.run_single_clone_model:
    run_type_shell_cmd = "python {params.script} -i {input} -o {output}"

else:
    run_type_shell_cmd = "mkdir -p {output}; touch {output}/full.txt; touch {output}/normal.txt;"


checkpoint write_run_type_sentinels:
    input:
        config.clone_cn_file,
    output:
        temp(directory(config.run_type_sentinel_dir)),
    params:
        script=workflow.source_path("scripts/write_run_type_sentinels.py"),
    conda:
        "envs/python.yaml"
    # noinspection SmkUnusedLogFile
    log:
        config.log_dir.joinpath("write_run_type_sentinels.log"),
    shell:
        "({}) >{{log}} 2>&1".format(run_type_shell_cmd)


rule run_cfclone:
    input:
        c=config.clone_cn_file,
        i=config.ctdna_file,
    output:
        f=config.fit_template,
        d=directory(config.exec_dir),
    params:
        c=config.num_chains,
        r=config.num_rounds,
        v=config.num_chains_vi,
        o=config.get_cfclone_outlier_args,
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
        "--clone-cnv-file {input.c} "
        "--in-file {input.i} "
        "--out-file {output.f} "
        "--exec-dir {output.d} "
        "--num-threads {threads} "
        "--num-chains {params.c} "
        "--num-chains-vi {params.v} "
        "--num-rounds {params.r} "
        "{params.o} "
        "{params.rt}) >{log} 2>&1"


rule write_ancestral_prevlances:
    input:
        i=str(config.fit_template).format(run_type="full"),
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
        str(config.fit_template).format(run_type="full"),
    output:
        config.dominance_prob_file,
    conda:
        "envs/cfclone.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.dominance_prob_file),
    shell:
        "(cfclone write-dominance-prob -i {input} -o {output}) >{log} 2>&1"


rule write_pairwise_ranks_file:
    input:
        str(config.fit_template).format(run_type="full"),
    output:
        config.pairwise_ranks_file,
    conda:
        "envs/cfclone.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.pairwise_ranks_file),
    shell:
        "(cfclone write-pairwise-ranks -i {input} -o {output}) >{log} 2>&1"


rule write_parameter_summaries_file:
    input:
        str(config.fit_template).format(run_type="full"),
    output:
        config.parameter_summaries_file,
    conda:
        "envs/cfclone.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.parameter_summaries_file),
    shell:
        "(cfclone write-parameter-summaries -i {input} -o {output}) >{log} 2>&1"


rule write_summary_file:
    input:
        str(config.fit_template).format(run_type="full"),
    output:
        config.summary_file,
    conda:
        "envs/cfclone.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.summary_file),
    shell:
        "(cfclone write-summary -i {input} -o {output}) >{log} 2>&1"


rule write_tumour_content_file:
    input:
        str(config.fit_template).format(run_type="full"),
    output:
        config.tumour_content_file,
    conda:
        "envs/cfclone.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.tumour_content_file),
    shell:
        "(cfclone write-tumour-content -i {input} -o {output}) >{log} 2>&1"


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


def get_runtypes(wildcards):
    checkpoint_output = checkpoints.write_run_type_sentinels.get(**wildcards).output[0]
    return glob_wildcards(os.path.join(checkpoint_output, "{run_type}.txt")).run_type


rule merge_evidence:
    input:
        expand(config.run_type_evidence_template, run_type=get_runtypes),
    output:
        config.evidence_file,
    params:
        script=workflow.source_path("scripts/merge_evidence.py"),
    conda:
        "envs/python.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.evidence_file),
    shell:
        "(python {params.script} -i {input} -o {output}) >{log} 2>&1"


rule plot_clone_prevalences:
    input:
        d=config.ancestral_prevalence_file,
        t=config.clone_prevalence_tree_file,
    output:
        config.clone_prevalences_plot,
    params:
        script=workflow.source_path("scripts/plot_clone_prevalences.py"),
    conda:
        "envs/plot.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.clone_prevalence_tree_file),
    shell:
        "(python {params.script} -d {input.d} -t {input.t} -o {output}) >{log} 2>&1"


rule plot_fit:
    input:
        config.parameter_summaries_file,
    output:
        config.fit_plot,
    params:
        script=workflow.source_path("scripts/plot_fit.py"),
    conda:
        "envs/plot.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.fit_plot),
    shell:
        "(python {params.script} -i {input} -o {output}) >{log} 2>&1"


rule plot_pairwsie_ranks:
    input:
        i=config.pairwise_ranks_file,
        t=config.clone_tree_file,
    output:
        config.pairwise_ranks_plot,
    params:
        script=workflow.source_path("scripts/plot_pairwise_ranks.py"),
    conda:
        "envs/plot.yaml"
    group:
        "post_process"
    log:
        config.get_log_file(config.pairwise_ranks_plot),
    shell:
        "(python {params.script} -i {input.i} -t {input.t} -o {output}) >{log} 2>&1"


rule save_run_configuration:
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
