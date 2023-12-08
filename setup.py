import setuptools

setuptools.setup(
    name="striplines",
    version="0.0.2",
    author="Will Henderson", 
    description="Calculates magnetic fields produced by striplines in real and fourier space", 
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3", 
        "License :: OSI Approved :: MIT License"
    ], 
    python_requires= '>3.0',
    py_modules=["striplines"],
    package_dir={'':'.'}
)