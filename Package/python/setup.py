from setuptools import setup,find_packages

setup(
        ### information about package and author.
        name = 'firstmatrixc',
        version = '0.2.1',
        author = 'shihua',
        author_email = "15021408795@163.com",
        python_requires = ">=3.9.13",
        license = "BSD 3 clause",

        ### source codes and dependencies
        packages = find_packages(),
        include_package_data = True,
        description = 'firstmatrixc is a ptyhon interface package that uses ctypes to encapsulate the C language matrix computation library.'
        # install_requires = ['ctypes','os'],
     
)