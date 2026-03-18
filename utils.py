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
        return self.config.get("num_chains", 8)

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
        return range(self.num_restarts)

    # Directories
    @property
    def out_dir(self):
        return Path(self.config["out_dir"])

    @property
    def pipeline_dir(self):
        return Path(self.config["pipeline_dir"])

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
    def restart_dir(self):
        return self.out_dir.joinpath("restart_{seed}")

    @property
    def ancestral_prevalence_file(self):
        return self.restart_dir.joinpath("tables", "ancestral_prevalence.tsv.gz")

    @property
    def dominance_prob_file(self):
        return self.restart_dir.joinpath("tables", "dominance_prob.tsv")

    @property
    def clone_prevalence_tree_file(self):
        return self.restart_dir.joinpath("trees", "prevalence_tree.json")

    @property
    def clone_prevalences_plot(self):
        return self.restart_dir.joinpath("plots", "clone_prevalences.pdf")

    @property
    def evidence_file(self):
        return self.restart_dir.joinpath("tables", "evidence.tsv")

    @property
    def exec_dir_template(self):
        return self.restart_dir.joinpath("fit", "{run_type}")

    @property
    def fit_template(self):
        return self.restart_dir.joinpath("fit", "{run_type}.h5")

    @property
    def fit_plot(self):
        return self.restart_dir.joinpath("plots", "fit.pdf")

    @property
    def pairwise_ranks_file(self):
        return self.restart_dir.joinpath("tables", "pairwise_ranks.tsv")

    @property
    def pairwise_ranks_plot(self):
        return self.restart_dir.joinpath("plots", "pairwise_ranks.pdf")

    @property
    def parameter_summaries_file(self):
        return self.restart_dir.joinpath("tables", "parameter_summaries.tsv.gz")

    @property
    def run_type_evidence_template(self):
        return self.tmp_dir.joinpath("evidence", "{seed}", "{run_type}.csv")

    @property
    def summary_file(self):
        return self.restart_dir.joinpath("tables", "summary.tsv")

    @property
    def tumour_content_file(self):
        return self.restart_dir.joinpath("tables", "tumour_content.tsv")

    @property
    def experiment_configuration(self):
        return self.out_dir.joinpath("config.yaml")

    @property
    def merged_evidence_file(self):
        return self.out_dir.joinpath("evidence.tsv")

    @property
    def merged_prevalence_file(self):
        return self.out_dir.joinpath("prevalence.tsv")

    @property
    def merged_summary_file(self):
        return self.out_dir.joinpath("summary.tsv")

    @property
    def merged_tumour_content_file(self):
        return self.out_dir.joinpath("tumour_content.tsv")

    @property
    def pipeline_files(self):
        result = []

        for s in self.seeds:
            result.append(self.merged_tumour_content_file)

            result.append(self.merged_evidence_file)

            result.append(self.merged_summary_file)

            result.append(str(self.fit_plot).format(seed=s))

            result.append(str(self.clone_prevalences_plot).format(seed=s))

            result.append(str(self.pairwise_ranks_plot).format(seed=s))

        return result

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
        if wildcards.run_type == "normal":
            return "--only-normal"
        else:
            return ""

    def _get_relative_path(self, template):
        try:
            rel_path = template.relative_to(self.pipeline_dir)
            parent = "working"
        except ValueError:
            rel_path = template.relative_to(self.out_dir)
            parent = "output"
        return parent, rel_path
