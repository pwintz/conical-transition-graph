#!/bin/env python
"""
CTG Example
"""

# import os
# import argparse
# import math
import numpy as np
# from polytope import polytope.Polytope as poly
import polytope as poly
from polytope import Polytope
# from poly import Polytope

print(dir(poly))

# Type hinting.
# from typing import List, Set, Dict, Tuple, Union, Callable

def main():
    """
    main function
    """
    # pylint: disable=locally-disabled, invalid-name
    A = np.array([[1.0, 0.0],
              [0.0, 1.0],
              [-1.0, -0.0],
              [-0.0, -1.0]])
    # pylint: enable=locally-disabled, invalid-name

    b = np.array([2.0, 1.0, 0.0, 0.0])

    p = Polytope(A, b)
    print(p)
    p.plot()


if __name__ == "__main__":
    main()
