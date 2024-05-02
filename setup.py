from setuptools import setup, find_packages

setup(
    name='artifacts_remove',  # Replace 'mypackage' with the name of your package
    version='0.1.0',  # Replace '0.1.0' with the current version of your package
    author='Yinzhao Wang',  # Replace 'Your Name' with your name
    author_email='wyz2016828@example.com',  # Replace with your email
    description='A package for detect and removal of artifacts in Spatial Transcriptomics Tissue samples',  # A short description
    long_description=open('README.md').read(),  # Long description read from the README.md
    long_description_content_type='text/markdown',  # This is important for formatting the long description
    url='https://github.com/KummerfeldLab/artifactDetection',  # Link to your package's GitHub repo
    packages=find_packages(),  # Finds all python modules in the directory automatically
    install_requires=[
        'numpy', 
        'pandas', 
        'copy'
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',  # Change as necessary
        'License :: OSI Approved :: MIT License',  # Choose the appropriate license
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    keywords='artifacts removal',  # Keywords for your project
    python_requires='>=3.6'  # Minimum version requirement of the Python runtime
)
