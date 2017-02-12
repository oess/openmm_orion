from __future__ import unicode_literals
"""
Copyright (C) 2016 OpenEye Scientific Software
"""
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube
from PlatformTestCubes.cubes import PlatformTestCube

job = WorkFloe("OpenMMPlatforms")

job.description = """
**Check available OpenMM Platforms**
Based on OpenMM SimTK installation check script
"""

job.classification = [
    ["OpenMM", "Platforms"]
]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ofs = OEMolOStreamCube("ofs")

# Promotes the parameter to something we can specify from the command line as "--ifs=..."
ifs.promote_parameter("data_in", promoted_name="ifs")

# this is hardwiring the filename to the molecules coming out of ofs
ofs.set_parameters(data_out="openmmPlatformCheck.oeb")

# the name of the object has to match the string: this is the name of myself
platformTester = PlatformTestCube("platformTester")

job.add_cubes(ifs, ofs, platformTester)
ifs.success.connect(platformTester.intake)
platformTester.success.connect(ofs.intake)


if __name__ == "__main__":
    job.run()
