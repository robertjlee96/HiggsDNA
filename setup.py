import setuptools
import yaml


def get_dependencies(env_yaml_file):
    """Scan a YAML environment file to get a list of dependencies"""
    with open(env_yaml_file, "r") as f:
        environment = yaml.safe_load(f)
    dependencies = []
    for dep in environment["dependencies"]:
        if not dep.startswith("python"):
            dependencies.append(dep)
    return dependencies


setuptools.setup(
    name="higgs_dna",
    packages=[
        "higgs_dna",
        "higgs_dna/systematics",
        "higgs_dna/tools",
        "higgs_dna/utils",
        "higgs_dna/workflows",
    ],
    scripts=["scripts/runner.py"],
    install_requires=get_dependencies("environment.yml"),
    extras_require={"dev": ["pytest>=3.7", "black>=19.10b0", "flake8>=3.8.4"]},
    # for metaconditions and JSON files in general
    include_package_data=True,
    python_requires=">=3.6",
)
