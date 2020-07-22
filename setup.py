import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="hyphalnet", # Replace with your own username
    version="0.0.1",
    author="Sara Gosline",
    author_email="sara.gosline@pnnl.gov",
    description="HyphalNet combines steiner tree reductions with multigraph community detection",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/sgosline/hyphalnet",
    packages=setuptools.find_packages(),
    install_requires=['numpy', 'pandas', 'matplotlib', 'goenrich', 'leidenalg', 'python_igraph'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
    package_data={'':['data/*txt','data/*obo','data/*gz']},
    python_requires='>=3.6')
