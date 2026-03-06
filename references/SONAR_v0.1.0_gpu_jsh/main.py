# import for system
import os
import sys
import time
import cupy as cp
from cupy.cuda import Device
# import for methods and functions
import initialization as init
import pre_process as pre
import computation as comp
import post_process as post
# import finalization as fin


def reset_gpu():
    with cp.cuda.Device(0):
        cp.get_default_pinned_memory_pool().free_all_blocks()
        cp._default_memory_pool.free_all_blocks()


def main(args, logs):
    # =========================================================================
    # 1. Pre-process
    # =========================================================================
    # Initialize and set variables.
    i1 = 0  # First index of the process
    i1 += 1
    logs.info('==== %i.    Pre-process' % i1)
    objs = pre.do_process(i1, args, logs)

    # =========================================================================
    # 2. Computation.
    # =========================================================================
    # Initialize and set variables.
    i1 += 1
    logs.info('==== %i.    Computation' % i1)
    comp.do_process(i1, args, logs, objs)

    # =========================================================================
    # 3. Post-process
    # =========================================================================
    i1 += 1
    logs.info('==== %i.    Post-process' % i1)
    post.do_process(objs)


# =============================================================================
#  Run main routine.
# =============================================================================
# this is only accessed if running from command prompt
if __name__ == '__main__':

    # GPU 메모리 초기화
    reset_gpu()

    # Parse arguments.
    args = init.parse_arguments()
    # Set loggers.
    logs = init.set_loggers(args)

    # Explain information of the code at the initial step.
    init.explain_information(logs)

    # Get wall-clock time at the initial step.
    # tic = time.time()

    # Run main routine without profiling procedure
    main(args, logs)

    # Get wall-clock time at the final step.
    # toc = time.time()

    # print(f"computing time = {toc - tic}")

    # Explain information of the code at the final step.
    # fin.explain_information(logs, toc-tic)

    # reset_gpu()
