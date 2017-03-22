
from io import StringIO

from floe.api import(
    OEMolComputeCube, parameter, BinaryOutputPort, SourceCube
)

from PlatformTestCubes import testInstallation
from PlatformTestCubes.benchmarking import run_platform_benchmarks


# For parallel, import and inherit from ParallelOEMolComputeCube
class PlatformTestCube(OEMolComputeCube):
    """
    Runs a copy of OpenMM's simtk.installation script

    Checks available OpenMM platforms
    """
    title = "OpenMM Platform Check"
    description = """
    *OpenMM Platform Check*
    Checks to see which OpenMM Platforms are available amongst CPU, OpenCL, and CUDA
    """
    classification = [
        ["OpenMM", "PlatformCheck"],
    ]
    tags = [tag for lists in classification for tag in lists]

    # the string in the first argument has to be the same as the name of the object
    # note that right now, the pdb file has to reside in the same directory as this file
    pdbFileName = parameter.StringParameter(
        "pdbFileName",
        title="pdb File Name",
        description="name of pdb file to use in OpenMM Platform Check",
        default='test.pdb',
    )

    def process(self, mol, port):
        output = testInstallation.run_tests(self.args.pdbFileName)
        self.log.info(output)
        with open('openmmPlatformCheck.txt', 'w') as ofs:
            ofs.write(output)
        self.emit(mol)


class BenchmarkCube(SourceCube):

    title = "OpenMM BenchmarkCube"

    description = """
        Cube that performs a benchmark of OpenMM on all of the different
        platforms that are available and outputs the byte string
        resulting from the benchmarks line by line to its success port
    """

    tags = [["OpenMM", "Benchmarking"]]

    success = BinaryOutputPort("success")

    cutoff = parameter.DecimalParameter("cutoff", default=0.9)
    seconds = parameter.IntegerParameter("seconds", default=60)
    polarization = parameter.StringParameter(
        "polarization",
        default="mutual",
        choices=["direct", "extrapolated", "mutual"]
    )
    amoeba_target_epsilon = parameter.DecimalParameter(
        "amoeba_target_epsilon",
        default=1e-5,
        title="Amoeba Mutual Induced Target Epsilon"
    )
    use_heavy_hydrogens = parameter.BooleanParameter(
        "use_heavy_hydrogens",
        default=False,
        title="Use Heavy Hydrogens"
    )
    precision = parameter.StringParameter(
        "precision",
        default="single",
        choices=["single", "mixed", "double"]
    )

    def __iter__(self):
        stream = StringIO()
        stream.write("Benchmarking Results:\n")
        run_platform_benchmarks(self.args, stream=stream)
        stream.flush()
        stream.seek(0)
        output = stream.readline()
        while len(output):
            self.log.info(output)
            yield output.encode("utf-8")
            output = stream.readline()
