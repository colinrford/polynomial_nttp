# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: 2025-2026 Colin Ford
#
# Conan 2.x recipe for lam.polynomial.nttp.
#
# Local development:
#   conan create . --profile <your-profile>
#
# Consumer projects can then declare:
#   requires = "lam_polynomial_nttp/<version>"   # version comes from VERSION
#
# This recipe is supplementary to the project's primary CMake+CPS distribution
# path — see CMakeLists.txt and the install(PACKAGE_INFO) block for the
# canonical install layout. The recipe delegates entirely to CMake.

from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, cmake_layout
from conan.tools.files import copy, load
import os

class LamPolynomialNttpConan(ConanFile):
    name = "lam_polynomial_nttp"
    license = "AGPL-3.0-or-later"
    author = "Colin Ford"
    url = "https://github.com/colinrford/polynomial_nttp"
    description = (
        "Compile-time polynomial algebra with non-type-template-parameter "
        "coefficients for the lam project (C++23 modules)."
    )
    topics = ("c++23", "modules", "polynomial", "nttp", "constexpr", "lam")

    settings = "os", "compiler", "build_type", "arch"
    package_type = "static-library"

    def set_version(self):
        # Single source of truth: the top-level VERSION file, which CMakeLists
        # also reads. Editing VERSION updates both the CMake project version and
        # this recipe — they never drift.
        self.version = load(
            self, os.path.join(self.recipe_folder, "VERSION")
        ).strip()

    def requirements(self):
        # lam.polynomial.nttp imports lam.concepts (identity traits, algebraic
        # concepts). Mirrors target_link_libraries(... lam_concepts::concepts).
        # Range, not exact pin: concepts uses date-based versions (0.1.YYMMDD)
        # that don't encode semver compatibility, so float on any 0.1.* at or
        # after this date and only re-pin at a 0.2 (minor) bump. Avoids editing
        # this recipe every time concepts bumps a date. A lockfile pins the
        # concrete version for reproducible builds.
        self.requires("lam_concepts/[>=0.1.260530 <0.2]")

    exports_sources = (
        "VERSION",
        "CMakeLists.txt",
        "polynomial_nttp_config.cppm.in",
        "src/*",
        "cmake/*",
        "LICENSE",
        "README.md",
    )

    def layout(self):
        cmake_layout(self)

    def validate(self):
        cppstd = self.settings.compiler.cppstd
        if cppstd is not None:
            std = int(str(cppstd).replace("gnu", ""))
            if std < 23:
                raise Exception(
                    "lam_polynomial_nttp requires C++23 (compiler.cppstd >= 23)."
                )

    def generate(self):
        tc = CMakeToolchain(self)
        # Mirror what the project's CMakeLists already assumes.
        tc.cache_variables["CMAKE_CXX_STANDARD"] = "23"
        tc.cache_variables["CMAKE_CXX_SCAN_FOR_MODULES"] = "ON"
        # tests/examples/benchmarks aren't in exports_sources — don't descend
        # into them during the package build.
        tc.cache_variables["LAM_POLYNOMIAL_NTTP_BUILD_TESTS"] = "OFF"

        # C++23 module deps: do NOT use a Conan deps generator here. Neither
        # CMakeDeps nor CMakeConfigDeps reproduces an imported target's
        # `FILE_SET CXX_MODULES`, so `import lam.concepts;` can't find the
        # module interface. Instead, point find_package at each dependency's
        # OWN installed config inside its Conan package folder — the real
        # lam_conceptsConfig.cmake carries the module file set. CMakeLists
        # already calls find_package(lam_concepts) via lam_find_dependency.
        prefix_paths = [
            self.dependencies[dep].package_folder
            for dep in ("lam_concepts",)
        ]
        tc.cache_variables["CMAKE_PREFIX_PATH"] = ";".join(prefix_paths)
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()
        copy(
            self,
            "LICENSE",
            src=self.source_folder,
            dst=os.path.join(self.package_folder, "licenses"),
        )

    def package_info(self):
        # Match the CMake export names so find_package(lam_polynomial_nttp) and
        # the Conan-generated CMakeDeps both yield the same target. Path A (per-
        # submodule primary CPS) namespaces this as
        # lam_polynomial_nttp::polynomial_nttp.
        self.cpp_info.set_property("cmake_file_name", "lam_polynomial_nttp")
        self.cpp_info.set_property(
            "cmake_target_name", "lam_polynomial_nttp::polynomial_nttp"
        )
        self.cpp_info.libs = ["lam_polynomial_nttp"]
        # CPS metadata lives under <pkg>/lib/cps for CMake 4.3+ consumers that
        # read CPS in addition to the legacy *Config.cmake.
        self.cpp_info.builddirs = [
            os.path.join("lib", "cmake", "lam_polynomial_nttp"),
            os.path.join("lib", "cps"),
        ]
