"""Навчальний стенд дискретного логарифма."""
from .core import BabyTable, Result, bsgs, brute_force, interval_bsgs, pohlig_hellman, pollard_rho
from .number_theory import crt, factor, is_prime, multiplicative_order, primitive_root

__all__ = ['BabyTable','Result','bsgs','brute_force','interval_bsgs','pohlig_hellman',
           'pollard_rho','crt','factor','is_prime','multiplicative_order','primitive_root']
