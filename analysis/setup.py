from setuptools import setup

setup(
    
    name = "pyouch.analysis",
    verison = "0.0.1",
    description = "tools for analyze results from ou model fitting",
    url = "https://github.com/Seongeun02/pyouch.git",
    author = "SeongeunYun",
    author_email = "lizsarah@snu.ac.kr",
    install_requires = [
        "numpy", "pandas", "matplotlib.pyplot", "os", "seaborn",
        "scipy", "sklearn", "statsmodels", "warnings", "re"
    ]

)
