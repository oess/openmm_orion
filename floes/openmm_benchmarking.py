from __future__ import unicode_literals
"""
Copyright (C) 2016 OpenEye Scientific Software
"""
from floe.api import WorkFloe, FileOutputCube
from PlatformTestCubes.cubes import BenchmarkCube

job = WorkFloe("OpenMM Benchmarking")

job.description = """
Performs Benchmarking upon all available Platforms
Based on OpenMM SimTK Benchmarking script
"""

job.classification = [
    ["OpenMM", "Platforms", "Benchmarking"]
]
job.tags = [tag for lists in job.classification for tag in lists]

benchmark_cube = BenchmarkCube("benchmark_cube")
ofs = FileOutputCube("ofs")

ofs.set_parameters(name="Orion_OpenMM_Benchmarks.txt")


job.add_cubes(benchmark_cube, ofs)
benchmark_cube.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
