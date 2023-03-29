import argparse
import json
import subprocess

def launch_tower(pbmc_run_params, barnyard_run_params, branch, repo_name):
    tower_download = "curl -L https://github.com/seqeralabs/tower-cli/releases/download/v0.7.3/tw-0.7.3-linux-x86_64 > tw"
    subprocess.run(tower_download, shell=True)
    
    make_tower_executable = "chmod +x ./tw"
    subprocess.run(make_tower_executable, shell=True)
    
    with open("pbmcRunParams.json", "w") as f:
        json.dump(pbmc_run_params, f)

    with open("barnyardRunParams.json", "w") as f:
        json.dump(barnyard_run_params, f)

    pipeline = "nf-rna" if repo_name == "nf-rna" else "ScaleRNA"
    
    tower_launch_pbmc = f"./tw launch {pipeline} --params-file=pbmcRunParams.json --revision {branch} --name=pbmc_test_pr{pbmc_run_params['pr_number']} --workspace=49873529197394"
    subprocess.run(tower_launch_pbmc, shell=True)

    tower_launch_barnyard = f"./tw launch {pipeline} --params-file=barnyardRunParams.json --revision {branch} --name=barnyard_test_pr{barnyard_run_params['pr_number']} --workspace=49873529197394"
    subprocess.run(tower_launch_barnyard, shell=True)

def build_pbmc_test_params(pr_number):
    return {"outDir": "s3://scale.bioinfo/testing/datasets/pbmc/test_out",
            "genome": "s3://scale.bioinfo/testing/datasets/pbmc/genome.json",
            "runFolder": "s3://scale.pub/testData/rna/PBMC/221104_NB551228_0182_AHCCFKBGXN",
            "samples": "s3://scale.bioinfo/testing/datasets/pbmc/samples.csv",
            "pr_number": pr_number}

def build_barnyard_test_params(pr_number):
    return {"outDir": "s3://scale.bioinfo/testing/datasets/barnyard/test_out",
            "genome": "s3://scale.bioinfo/genomes/grch38_mm39/grch38_mm39.json",
            "fastqDir": "s3://scale.bioinfo/testing/datasets/barnyard/fastq_files",
            "samples": "s3://scale.bioinfo/testing/datasets/barnyard/samples.csv",
            "internalReport": True,
            "pr_number": pr_number}

def main():
    parser = argparse.ArgumentParser(description="Launch test pipeline on tower")
    parser.add_argument("--branch", type=str, help="Git branch to run pipeline on")
    parser.add_argument("--pr_number", type=str, help="Git pull request number")
    parser.add_argument("--repo_name", type=str, help="Git repo name")

    args = parser.parse_args()

    pbmc_run_params = build_pbmc_test_params(args.pr_number)

    barnyard_run_params = build_barnyard_test_params(args.pr_number)

    launch_tower(pbmc_run_params, barnyard_run_params, args.branch, args.repo_name)

if __name__ == '__main__':
    main()
