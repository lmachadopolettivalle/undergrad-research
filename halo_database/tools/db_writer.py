#!/usr/bin/env python2.7
import matplotlib
import sys
from tangos import parallel_tasks, core
from tangos.tools.property_writer import PropertyWriter


def run_dbwriter(argv):
    parallel_tasks.mpi_sync_db(core.get_default_session())
    writer = PropertyWriter()
    writer.parse_command_line(argv)
    writer.run_calculation_loop()


if __name__ == "__main__":
    #matplotlib.use('agg')
    parallel_tasks.launch(run_dbwriter, 2, [sys.argv[1:]])

