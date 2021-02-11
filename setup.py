import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="weave-weave_ifu_workflow",
    version="0.0.1",
    description="Workflow for configuring weave fields",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages = setuptools.find_packages(),
    python_requires='>=3',
    install_requires=['numpy', 'astropy', 'astropy-healpix', 'matplotlib']
)