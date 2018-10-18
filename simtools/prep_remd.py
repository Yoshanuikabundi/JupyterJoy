import sys
import os
import mdtraj as md
from math import ceil, sqrt
from enum import Enum
from typing import List, Optional


def calc_temps_optimal(
    min_temp: float,
    num_reps: int,
    const: float,
    scalefirst: float = 1.0
) -> List[float]:
    """Compute a temperature ladder with num_reps rungs starting at min_temp

    Method from Prakash, Barducci and Parrinello, 2011, JCTC 7 (7) 2025-2027
    https://pubs.acs.org/doi/full/10.1021/ct200208h
    Const is c/a in eq 5
    The size of the first step on the ladder will be scaled by scalefirst."""
    temps = [min_temp]
    for i in range(1, num_reps):
        t_im1 = temps[-1]
        b_im1 = 1 / t_im1
        b_i = (b_im1 - sqrt(const * b_im1))
        t_i = 1 / b_i
        if i == 1:
            t_i = (t_i - t_im1) * scalefirst + t_im1
        temps.append(t_i)
    return temps


def calc_optimal_const(
    min_temp: float,
    max_temp: float,
    num_reps: int,
    lbound: float = 0.0,
    ubound: Optional[float] = None,
    tolerance: Optional[float] = None,
    scalefirst: float = 1.0
) -> float:
    """Find an appropriate constant for calc_temps_optimal()

    Binary search between lbound and ubound for a constant that reproduces
    the desired min/max/num within tolerance. Overkill but easy search
    problem so who cares. Defaults for lbound, ubound and tolerance
    should return a constant that gives min/max/num to within machine
    precision for almost any combination of min/max/num. The size fo the
    first step on the ladder will be scaled by scalefirst."""
    if tolerance is None:
        tolerance = sys.float_info.epsilon * max_temp

    if ubound is None:
        ubound = 1.0 / (min_temp * num_reps ** 2)  # This seems to work,
        # don't know why. There's some room to wiggle the 3.4 but it's
        # within unit of optimal, eg 4.0 fails for every combo I tried and
        # 3.0 is probably too low. Ideally we'd choose a bound such that
        # calc_temps_optimal(min_temp, num_reps, const)[-1] == min(
        #     sys.float_info.max,
        #     1.0 / sys.float_info.min
        # ) # Assuming min() were computed with infinite precision
        # I think that would need it's own binary search, which is totally
        # doable, but that's overkill.

    if min_temp >= max_temp:
        raise ValueError('max_temp must be strictly greater than min_temp')

    # Make sure the specified bounds are OK
    if calc_temps_optimal(min_temp, num_reps, lbound)[-1] > max_temp:
        raise ValueError(
            'Lower bound on const implies maximum temperature higher than max!'
        )
    if calc_temps_optimal(min_temp, num_reps, ubound)[-1] < max_temp:
        raise ValueError(
            'Upper bound on const implies maximum temperature lower than max!'
        )

    while True:
        # Try a const midway between the bounds
        const = (ubound + lbound) / 2.0
        max_temp_const = calc_temps_optimal(
            min_temp,
            num_reps,
            const,
            scalefirst=scalefirst
        )[-1]

        # If it's within tolerance, return const
        if abs(max_temp_const - max_temp) < tolerance:
            return const
        # Otherwise, update the appropriate bound and try again
        elif max_temp_const > max_temp:
            ubound = const
        elif max_temp_const < max_temp:
            lbound = const
        else:
            raise TypeError('const value {const} makes no sense (probs NAN)')


def calc_temps(
    min_temp: float,
    max_temp: float,
    num_reps: int,
    method: str = 'OPTIMAL',
    **kwargs
) -> List[float]:
    """Calculate a temperature ladder"""
    if method == 'OPTIMAL':
        const = calc_optimal_const(min_temp, max_temp, num_reps, **kwargs)
        return calc_temps_optimal(min_temp, num_reps, const, **kwargs)
    elif kwargs:
        raise ValueError('Extra keyword arguments are only for method OPTIMAL')
    elif method == 'GEOMETRIC':
        ratio = (max_temp / min_temp) ** (1 / (num_reps - 1))
        return [min_temp * ratio ** i for i in range(num_reps)]
    elif method == 'LINEAR':
        step = (max_temp - min_temp) / (num_reps - 1)
        return [min_temp + step * i for i in range(num_reps)]
    else:
        raise ValueError(f'Unrecognised method: {method}')
