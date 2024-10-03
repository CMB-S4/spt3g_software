import logging
import os
import re
import subprocess
import sys
import sysconfig
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
from setuptools.command.install_scripts import install_scripts
from setuptools.command.sdist import sdist


# A CMakeExtension does not need a source list
class CMakeExtension(Extension):
    def __init__(self, name):
        super().__init__(name, sources=[])


class CMakeBuildExt(build_ext):
    def initialize_options(self):
        super().initialize_options()
        self.source_dir = Path(".").resolve()
        self.build_dir = Path(os.getenv("BUILD_DIR", "./build")).resolve()

        cmake_args = []

        # Adding CMake arguments set as environment variable
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        # sensible defaults
        if not any(["Python_ROOT_DIR" in a for a in cmake_args]):
            pyroot = sysconfig.get_config_var("prefix")
            cmake_args += [f"-DPython_ROOT_DIR={pyroot}"]
        if not any(["Python_EXECUTABLE" in a for a in cmake_args]):
            cmake_args += [f"-DPython_EXECUTABLE={sys.executable}"]

        # install libraries and headers into virtual environment too
        if "VIRTUAL_ENV" in os.environ:
            if not any(["CMAKE_INSTALL_PREFIX" in a for a in cmake_args]):
                envroot = os.getenv("VIRTUAL_ENV")
                cmake_args += [f"-DCMAKE_INSTALL_PREFIX={envroot}"]

        # pass version to C++ code
        if not any(["SPT3G_VERSION" in a for a in cmake_args]):
            cmake_args += [f"-DPIP_SPT3G_VERSION_FILE=ON"]
            version = self.distribution.metadata.version
            if re.match(r"^[0-9]+\.[0-9]+$", version):
                cmake_args += [f"-DSPT3G_VERSION={repr(version)}"]

        self.cmake_args = cmake_args
        self.cmake_bin = os.getenv("CMAKE_BIN", "cmake")
        self.cmake_done = False

    def build_extension(self, ext):
        if not self.cmake_done:
            # build once
            self.announce(f"Building library in {self.build_dir}", logging.INFO)
            if not self.build_dir.exists():
                self.build_dir.mkdir(parents=True, exist_ok=True)
            subprocess.run(
                [self.cmake_bin, self.source_dir, *self.cmake_args], cwd=self.build_dir, check=True
            )
            subprocess.run([self.cmake_bin, "--build", "."], cwd=self.build_dir, check=True)

            # update package directory
            if self.editable_mode:
                pkgdir = self.build_dir.relative_to(self.source_dir) / "spt3g"
                self.distribution.package_dir["spt3g"] = pkgdir
                build_py = self.get_finalized_command("build_py")
                build_py.package_dir["spt3g"] = pkgdir

            # update version file
            self.copy_file(self.source_dir / "cmake/package/version.py", self.build_dir)

            # trigger script installer for all python scripts
            self.distribution.scripts = [
                str(f) for f in self.build_dir.glob("bin/*") if f.is_symlink()
            ]
            self.cmake_done = True

        # add modules
        spt3g_lib = Path(self.build_lib) / "spt3g"
        if not spt3g_lib.exists():
            spt3g_lib.mkdir(parents=True)

        libname = ext.name.split(".")[-1]
        libfiles = [
            self.build_dir / "spt3g" / self.get_ext_filename(libname),
            self.build_dir / "spt3g" / f"{libname}.so",
        ]
        for f in libfiles:
            if f.exists():
                self.copy_file(f, spt3g_lib)
                break
        else:
            raise RuntimeError(f"Missing extension module {ext.name}")


class CMakeInstallScripts(install_scripts):
    def run(self):
        super().run()

        self.announce("Installing cmake programs", logging.INFO)
        install_dir = Path(self.install_dir)
        if not install_dir.exists():
            install_dir.mkdir(parents=True)

        for src in self.get_finalized_command("build_ext").build_dir.glob("bin/*"):
            self.copy_file(src, install_dir)


class CMakeInstall(install):
    def run(self):
        super().run()

        if "CIBUILDWHEEL" in os.environ:
            return

        self.announce("Installing cmake libraries and headers", logging.INFO)
        build_dir = self.get_finalized_command("build_ext").build_dir
        subprocess.run(["make", "install"], cwd=build_dir, check=True)


class ArchiveDist(sdist):
    def run(self):

        if Path(".git").exists():
            ar = subprocess.run(
                ["git", "archive", "HEAD", ".git_archival.txt"],
                check=True,
                capture_output=True,
            )
            subprocess.run(["tar", "-x"], input=ar.stdout, check=True)

        super().run()

        if Path(".git").exists():
            subprocess.run(["git", "checkout", ".git_archival.txt"], check=True)


# gather libraries
clibs = []
pdirs = {"spt3g": "cmake/package"}

for d in sorted(Path("./").glob("*/CMakeLists.txt")):
    lib = d.parent.name
    if (d.parent / "src").exists():
        clibs.append(f"spt3g._lib{lib}")
    if (d.parent / "python").exists():
        pdirs[f"spt3g.{lib}"] = d.parent / "python"
    elif (d.parent / "__init__.py").exists():
        pdirs[f"spt3g.{lib}"] = d.parent

setup(
    ext_modules=[CMakeExtension(lib) for lib in clibs],
    cmdclass={
        "build_ext": CMakeBuildExt,
        "install_scripts": CMakeInstallScripts,
        "install": CMakeInstall,
        "sdist": ArchiveDist,
    },
    packages=list(pdirs),
    package_dir=pdirs,
    exclude_package_data={k: ["CMakeLists.txt"] for k in pdirs},
)
