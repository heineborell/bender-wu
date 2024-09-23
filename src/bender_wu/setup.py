from setuptools import setup

requirements = ["sympy", "numpy", "symengine"]


requirements_dev = ["black", "flake8", "isort", "pre-commit"]


setup(
    name="bender_wu",
    version="0.1.0",
    description="computing QNMs using bender-wu algorithm",
    author="deniz",
    url="https://github.com/heineborell/bender-wu",
    packages=["bender_wu"],
    package_dir={"": "src"},
    install_requires=requirements,
    extra_require={"dev": requirements_dev},
)
