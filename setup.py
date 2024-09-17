import os
import subprocess
import sys
import sysconfig
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install


# A CMakeExtension needs a sourcedir instead of a file list.
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        super().__init__(name, sources=[])
        self.sourcedir = os.fspath(Path(sourcedir).resolve())


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        cmake_args = ["-DBUILD_WHEEL=yes"]

        # Adding CMake arguments set as environment variable
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        # sensible defaults
        if not any(["Python_ROOT_DIR" in a for a in cmake_args]):
            pyroot = sysconfig.get_config_var("prefix")
            cmake_args += [f"-DPython_ROOT_DIR={pyroot}"]

        build_temp = Path("wheel/build")
        if not build_temp.exists():
            build_temp.mkdir(parents=True)

        libfile = build_temp / "spt3g" / self.get_ext_filename(ext.name.split(".")[-1])
        if not libfile.exists():
            # build once
            subprocess.run(
                ["cmake", ext.sourcedir, *cmake_args], cwd=build_temp, check=True
            )
            subprocess.run(["cmake", "--build", "."], cwd=build_temp, check=True)

        # add modules
        spt3g_lib = Path(self.build_lib) / "spt3g"
        if not spt3g_lib.exists():
            spt3g_lib.mkdir(parents=True)
        self.copy_file(libfile, spt3g_lib)


# gather libraries
clibs = []
pdirs = {"spt3g": "./wheel/spt3g"}

for d in sorted(Path("./").glob("*/CMakeLists.txt")):
    lib = d.parent.name
    if (d.parent / "src").exists():
        clibs.append(lib)
    if (d.parent / "python").exists():
        pdirs[f"spt3g.{lib}"] = d.parent / "python"
    elif (d.parent / "__init__.py").exists():
        pdirs[f"spt3g.{lib}"] = d.parent

setup(
    ext_modules=[CMakeExtension(f"spt3g._lib{lib}") for lib in clibs],
    cmdclass={"build_ext": CMakeBuild},
    packages=list(pdirs),
    package_dir=pdirs,
    exclude_package_data={k: ["CMakeLists.txt"] for k in pdirs},
)
