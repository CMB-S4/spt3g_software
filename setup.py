import json
import os
import subprocess
import sys
import sysconfig
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
from setuptools.command.install_scripts import install_scripts


# A CMakeExtension needs a sourcedir instead of a file list.
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        super().__init__(name, sources=[])
        self.sourcedir = os.fspath(Path(sourcedir).resolve())


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        cmake_args = []

        # Adding CMake arguments set as environment variable
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        if not any(["BUILD_WHEEL" in a for a in cmake_args]):
            if "CI_BUILD" in os.environ:
                cmake_args += ["-DBUILD_WHEEL=yes"]

        # sensible defaults
        if not any(["Python_ROOT_DIR" in a for a in cmake_args]):
            pyroot = sysconfig.get_config_var("prefix")
            cmake_args += [f"-DPython_ROOT_DIR={pyroot}"]
        if not any(["Python_EXECUTABLE" in a for a in cmake_args]):
            cmake_args += [f"-DPython_EXECUTABLE={sys.executable}"]
        broot = Path("wheel/deps").resolve()
        if broot.exists():
            cmake_args += [f"-DBOOST_ROOT={broot}"]

        wbdir = Path(ext.sourcedir) / "wheel/build"
        if self.editable_mode or "CI_BUILD" in os.environ:
            build_temp = wbdir
        else:
            build_temp = Path(self.build_temp)
        if not build_temp.exists():
            build_temp.mkdir(parents=True)
        self.build_dir = build_temp

        libname = ext.name.split(".")[-1]
        libfiles = [
            build_temp / "spt3g" / self.get_ext_filename(libname),
            build_temp / "spt3g" / f"{libname}.so",
        ]
        if not any([f.exists() for f in libfiles]):
            # build once
            self.announce(f"Building library in {build_temp}")
            subprocess.run(
                ["cmake", ext.sourcedir, *cmake_args], cwd=build_temp, check=True
            )
            subprocess.run(["cmake", "--build", "."], cwd=build_temp, check=True)

        # add modules
        spt3g_lib = Path(self.build_lib) / "spt3g"
        if not spt3g_lib.exists():
            spt3g_lib.mkdir(parents=True)
        for f in libfiles:
            if f.exists():
                self.copy_file(f, spt3g_lib)
                break
        else:
            raise RuntimeError(f"Missing extension module {libname}")


class CMakeInstallScripts(install_scripts):
    def run(self):
        super().run()

        self.announce("Installing scripts")
        install_dir = Path(self.install_dir)
        if not install_dir.exists():
            install_dir.mkdir(parents=True)

        for src in self.get_finalized_command("build_ext").build_dir.glob("bin/*"):
            self.copy_file(src, install_dir)


class CMakeInstall(install):
    def run(self):
        super().run()

        if "CI_BUILD" in os.environ:
            return

        self.announce("Installing libraries and headers")
        build_dir = self.get_finalized_command("build_ext").build_dir
        subprocess.run(["make", "install"], cwd=build_dir, check=True)


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
    cmdclass={
        "build_ext": CMakeBuild,
        "install_scripts": CMakeInstallScripts,
        "install": CMakeInstall,
    },
    packages=list(pdirs),
    package_dir=pdirs,
    exclude_package_data={k: ["CMakeLists.txt"] for k in pdirs},
)
