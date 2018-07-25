import os
from glob import glob
from collection import (ImageCollection as GMOSPipeline, symbolic_link)

program_ids = ("GS-2017A-Q-66", "GS-2016A-Q-76", "GS-2016B-Q-80", )

linked_files = ("run.py", "utils.py", "Atlas.Arcturus.372_926nm.txt")

for program_id in program_ids:

    commands = []

    raw_data_dir = "../data/gemini/{}/raw".format(program_id)
    reduced_data_dir = "../data/gemini/{}/reduced_v2".format(program_id)

    # Download the data!

    # Create pipeline.
    pipeline = GMOSPipeline(glob(os.path.join(raw_data_dir, "*.fits")))

    # Identify associations.
    associations = pipeline.prepare(reduced_data_dir)

    # Begin reductions.
    for wd, association in associations.items():
        association.write(os.path.join(wd, "summary.fits"), overwrite=True)

        for path in linked_files:
            symbolic_link(path, wd)
    
        commands.extend([
            "cd {}".format(os.path.abspath(wd)),
            "mkiraf -f",
            "python run.py"
        ])


    commands.append("cd {}".format(os.path.dirname(__file__)))
    with open("{}.sh".format(program_id), "w") as fp:
        fp.write("\n".join(commands))


with open("reduce_all.sh", "w") as fp:
    fp.write("\n".join(["sh {}.sh".format(pid) for pid in program_ids]))