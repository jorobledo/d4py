from numpy.distutils.core import setup, Extension


setup(
        package_dir = {'d4py': 'd4py'},
        ext_modules=[
            Extension(name="pymc", sources=["d4py/mc.f"], extra_compile_args=["-c"])
            ]
        )
