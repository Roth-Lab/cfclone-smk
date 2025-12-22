import pandas as pd
import yaml

from pathlib import Path


class ConfigManager(object):
    def __init__(self, config):
        self.config = config

    @property
    def config_yaml(self):
        return yaml.dump(self.config)

    # Params
    @property
    def num_chains(self):
        return self.config.get("num_chains", 16)

    @property
    def num_chains_vi(self):
        return self.config.get("num_chains_vi", 8)

    @property
    def num_restarts(self):
        return self.config.get("num_restarts", 1)

    @property
    def num_rounds(self):
        return self.config.get("num_rounds", 10)

    @property
    def num_threads(self):
        return self.config.get("num_threads", 1)

    @property
    def seeds(self):
        return list(range(self.num_restarts))

    @property
    def run_single_clone_model(self):
        return self.config.get("run_single_clone_model", False)

    # Directories
    @property
    def out_dir(self):
        return Path(self.config.get("out_dir", "<out_dir>"))

    @property
    def restart_out_dir(self):
        return self.out_dir.joinpath("restart_{seed}")

    @property
    def pipeline_dir(self):
        return Path(self.config.get("pipeline_dir", "<pipeline_dir>"))

    @property
    def benchmark_dir(self):
        return self.pipeline_dir.joinpath("benchmark")

    @property
    def log_dir(self):
        return self.pipeline_dir.joinpath("log")

    @property
    def tmp_dir(self):
        return self.pipeline_dir.joinpath("tmp")

    # Input files
    @property
    def clone_cn_file(self):
        return Path(self.config["clone_cn_file"]).resolve()

    @property
    def clone_tree_file(self):
        return Path(self.config["clone_tree_newick"]).resolve()

    @property
    def ctdna_file(self):
        return Path(self.config["ctdna_file"]).resolve()

    # Pipeline files
    @property
    def ancestral_prevalence_file(self):
        return self.restart_out_dir.joinpath("tables", "ancestral_prevalence.tsv.gz")

    @property
    def dominance_prob_file(self):
        return self.restart_out_dir.joinpath("tables", "dominance_prob.tsv")

    @property
    def clone_prevalence_tree_file(self):
        return self.restart_out_dir.joinpath("trees", "prevalence_tree.json")

    @property
    def clone_prevalences_plot(self):
        return self.restart_out_dir.joinpath("plots", "clone_prevalences.pdf")

    @property
    def run_type_sentinel_dir(self):
        return self.tmp_dir.joinpath("run_type_sentinels")

    @property
    def evidence_file(self):
        return self.restart_out_dir.joinpath("tables", "evidence.tsv")

    @property
    def experiment_configuration(self):
        return self.out_dir.joinpath("config.yaml")

    @property
    def fit_template(self):
        return self.restart_out_dir.joinpath("fit", "{run_type}.h5")

    @property
    def fit_exec_dir_template(self):
        return self.restart_out_dir.joinpath("fit", "{run_type}")

    @property
    def fit_plot(self):
        return self.restart_out_dir.joinpath("plots", "fit.pdf")

    @property
    def merged_evidence_file(self):
        return self.out_dir.joinpath("merged_evidence.tsv")

    @property
    def merged_summary_file(self):
        return self.out_dir.joinpath("merged_summary.tsv")

    @property
    def merged_tumour_content_file(self):
        return self.out_dir.joinpath("merged_tumour_content.tsv")

    @property
    def pairwise_ranks_file(self):
        return self.restart_out_dir.joinpath("tables", "pairwise_ranks.tsv")

    @property
    def pairwise_ranks_plot(self):
        return self.restart_out_dir.joinpath("plots", "pairwise_ranks.pdf")

    @property
    def parameter_summaries_file(self):
        return self.restart_out_dir.joinpath("tables", "parameter_summaries.tsv.gz")

    @property
    def summary_file(self):
        return self.restart_out_dir.joinpath("tables", "summary.tsv")

    @property
    def tumour_content_file(self):
        return self.restart_out_dir.joinpath("tables", "tumour_content.tsv")

    @property
    def run_type_evidence_template(self):
        return self.tmp_dir.joinpath("restart_{seed}", "evidence", "{run_type}.csv")

    @property
    def pipeline_files(self):
        def format_file_name(x):
            return str(x).format(seed=s)

        files = []

        for s in self.seeds:
            files.append(format_file_name(self.ancestral_prevalence_file))

            files.append(format_file_name(self.clone_prevalences_plot))

            files.append(format_file_name(self.fit_plot))

            files.append(format_file_name(self.pairwise_ranks_plot))

        files.append(self.experiment_configuration)

        files.append(self.merged_evidence_file)

        files.append(self.merged_summary_file)

        files.append(self.merged_tumour_content_file)

        return files

    # Helper functions
    def get_benchmark_file(self, template):
        parent, rel_path = self._get_relative_path(template)
        rel_path = rel_path.with_suffix(".log")
        return self.benchmark_dir.joinpath(parent, rel_path)

    def get_log_file(self, template):
        parent, rel_path = self._get_relative_path(template)
        rel_path = rel_path.with_suffix(".log")
        return self.log_dir.joinpath(parent, rel_path)

    @staticmethod
    def get_cfclone_run_type_args(wildcards):
        run_type = wildcards.run_type
        if run_type == "full":
            return ""
        elif run_type == "normal":
            return "--only-normal"
        else:
            clone = run_type.split("_")[-1]
            return "--use-clone {}".format(clone)

    def _get_relative_path(self, template):
        try:
            rel_path = template.relative_to(self.pipeline_dir)
            parent = "working"
        except ValueError:
            rel_path = template.relative_to(self.out_dir)
            parent = "output"
        return parent, rel_path
