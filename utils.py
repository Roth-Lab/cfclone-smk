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
    def num_rounds(self):
        return self.config.get("num_rounds", 10)

    @property
    def num_threads(self):
        return self.config.get("num_threads", 1)

    @property
    def run_single_clone_model(self):
        return self.config.get("run_single_clone_model", False)

    # Directories
    @property
    def out_dir(self):
        return Path(self.config.get("out_dir", "<out_dir>"))

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
        return self.out_dir.joinpath("tables", "ancestral_prevalence.tsv.gz")

    @property
    def dominance_prob_file(self):
        return self.out_dir.joinpath("tables", "dominance_prob.tsv")

    @property
    def clone_prevalence_tree_file(self):
        return self.out_dir.joinpath("trees", "prevalence_tree.json")

    @property
    def clone_prevalences_plot(self):
        return self.out_dir.joinpath("plots", "clone_prevalences.pdf")

    @property
    def run_type_sentinel_dir(self):
        return self.tmp_dir.joinpath("run_type_sentinels")

    @property
    def evidence_file(self):
        return self.out_dir.joinpath("tables", "evidence.tsv")

    @property
    def experiment_configuration(self):
        return self.out_dir.joinpath("config.yaml")

    @property
    def fit_template(self):
        return self.out_dir.joinpath("fit", "{run_type}.h5")
    
    @property
    def exec_dir(self):
        return self.out_dir.joinpath("fit", "{run_type}")

    @property
    def fit_plot(self):
        return self.out_dir.joinpath("plots", "fit.pdf")

    @property
    def pairwise_ranks_file(self):
        return self.out_dir.joinpath("tables", "pairwise_ranks.tsv")

    @property
    def pairwise_ranks_plot(self):
        return self.out_dir.joinpath("plots", "pairwise_ranks.pdf")

    @property
    def parameter_summaries_file(self):
        return self.out_dir.joinpath("tables", "parameter_summaries.tsv.gz")

    @property
    def summary_file(self):
        return self.out_dir.joinpath("tables", "summary.tsv")

    @property
    def tumour_content_file(self):
        return self.out_dir.joinpath("tables", "tumour_content.tsv")

    @property
    def run_type_evidence_template(self):
        return self.tmp_dir.joinpath("evidence", "{run_type}.csv")

    @property
    def pipeline_files(self):
        return [
            str(self.fit_template).format(run_type="full"),
            str(self.exec_dir).format(run_type="full"),
        ]

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
