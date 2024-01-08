# This file is part of tblite.
# SPDX-Identifier: LGPL-3.0-or-later
#
# tblite is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tblite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with tblite.  If not, see <https://www.gnu.org/licenses/>.
"""
Thin wrapper around CFFI extension module tblite.

This module mainly acts as a guard for importing the libtblite extension and
also provides some FFI based wappers for memory handling.
"""

import functools
import sys
from typing import Callable

import numpy as np

try:
    from ._libtblite import ffi, lib
except ImportError:
    raise ImportError("tblite C extension unimportable, cannot use C-API")

from .exceptions import TBLiteRuntimeError, TBLiteTypeError


def get_version() -> tuple:
    """Return the current API version from tblite.
    For easy usage in C the API version is provided as

    10000 * major + 100 * minor + patch

    For Python we want something that looks like a semantic version again.
    """
    version = lib.tblite_get_version()
    return (
        version // 10000,
        version % 10000 // 100,
        version % 100,
    )


def _delete_error(error) -> None:
    """Delete a tblite error handler object"""
    ptr = ffi.new("tblite_error *")
    ptr[0] = error
    lib.tblite_delete_error(ptr)


def new_error():
    """Create new tblite error handler object"""
    return ffi.gc(lib.tblite_new_error(), _delete_error)


def error_check(func):
    """Handle errors for library functions that require an error handle"""

    @functools.wraps(func)
    def handle_error(*args, **kwargs):
        """Run function and than compare context"""
        _err = new_error()
        value = func(_err, *args, **kwargs)
        if lib.tblite_check_error(_err):
            _message = ffi.new("char[]", 512)
            lib.tblite_get_error(_err, _message, ffi.NULL)
            raise TBLiteRuntimeError(ffi.string(_message).decode())
        return value

    return handle_error


@ffi.def_extern()
def logger_callback(message, nchar, data):
    """Custom logger callback to write output in a Python friendly way"""

    callback = ffi.from_handle(data)
    callback(ffi.unpack(message, nchar).decode())


def context_check(func):
    """Handle errors for library functions that require a context handle"""

    @functools.wraps(func)
    def handle_context_error(ctx, *args, **kwargs):
        """Run function and than compare context"""
        if isinstance(ctx, tuple):
            ctx, _ = ctx
        value = func(ctx, *args, **kwargs)
        if lib.tblite_check_context(ctx):
            _message = ffi.new("char[]", 512)
            lib.tblite_get_context_error(ctx, _message, ffi.NULL)
            raise TBLiteRuntimeError(ffi.string(_message).decode())
        return value

    return handle_context_error


def _delete_context(context) -> None:
    """Delete a tblite context handler object"""
    ptr = ffi.new("tblite_context *")
    ptr[0] = context
    lib.tblite_delete_context(ptr)


def new_context(color: bool = True, logger: Callable[[str], None] = print):
    """Create new tblite context handler object"""
    ctx = ffi.gc(lib.tblite_new_context(), _delete_context)
    handle = ffi.new_handle(logger)
    context_check(lib.tblite_set_context_logger)(ctx, lib.logger_callback, handle)
    if color and sys.stdout.isatty():
        context_check(lib.tblite_set_context_color)(ctx, 1)
    return ctx, handle


def _delete_structure(mol) -> None:
    """Delete molecular structure data"""
    ptr = ffi.new("tblite_structure *")
    ptr[0] = mol
    lib.tblite_delete_structure(ptr)


def new_structure(natoms, numbers, positions, charge, uhf, lattice, periodic):
    """Create new molecular structure data"""
    return ffi.gc(
        error_check(lib.tblite_new_structure)(
            natoms,
            numbers,
            positions,
            charge,
            uhf,
            lattice,
            periodic,
        ),
        _delete_structure,
    )


update_structure_geometry = error_check(lib.tblite_update_structure_geometry)
update_structure_charge = error_check(lib.tblite_update_structure_charge)
update_structure_uhf = error_check(lib.tblite_update_structure_uhf)


def _delete_table(table):
    """Delete a tblite data table object"""
    ptr = ffi.new("tblite_table *")
    ptr[0] = table
    lib.tblite_delete_table(ptr)


def new_table():
    """Create a tblite data table object"""
    return ffi.gc(lib.tblite_new_table(ffi.NULL), _delete_table)


table_set_double = error_check(lib.tblite_table_set_double)
table_set_int64_t = error_check(lib.tblite_table_set_int64_t)
table_set_bool = error_check(lib.tblite_table_set_bool)
table_set_char = error_check(lib.tblite_table_set_char)
table_add_table = error_check(lib.tblite_table_add_table)


def _delete_param(param):
    """Delete a tblite parametrization record object"""
    ptr = ffi.new("tblite_param *")
    ptr[0] = param
    lib.tblite_delete_param(ptr)


def new_param():
    """Create a tblite data table object"""
    return ffi.gc(lib.tblite_new_param(), _delete_param)


load_param = error_check(lib.tblite_load_param)
dump_param = error_check(lib.tblite_dump_param)
export_gfn2_param = error_check(lib.tblite_export_gfn2_param)
export_gfn1_param = error_check(lib.tblite_export_gfn1_param)
export_ipea1_param = error_check(lib.tblite_export_ipea1_param)


def _delete_result(result) -> None:
    """Delete a tblite result container object"""
    ptr = ffi.new("tblite_result *")
    ptr[0] = result
    lib.tblite_delete_result(ptr)


def new_result():
    """Create new tblite result container object"""
    return ffi.gc(lib.tblite_new_result(), _delete_result)


def copy_result(res):
    """Create new tblite result container object as copy of an existing result"""
    return ffi.gc(lib.tblite_copy_result(res), _delete_result)


def get_number_of_atoms(res) -> int:
    """Retrieve number of atoms from result container"""
    _natoms = ffi.new("int *")
    error_check(lib.tblite_get_result_number_of_atoms)(res, _natoms)
    return _natoms[0]


def get_number_of_orbitals(res) -> int:
    """Retrieve number of orbitals from result container"""
    _norb = ffi.new("int *")
    error_check(lib.tblite_get_result_number_of_orbitals)(res, _norb)
    return _norb[0]


def get_energy(res) -> float:
    """Retrieve energy from result container"""
    _energy = np.array(0.0)
    error_check(lib.tblite_get_result_energy)(
        res, ffi.cast("double*", _energy.ctypes.data)
    )
    return _energy


def get_energies(res):
    """Retrieve atom-resolved energies from result container"""
    _energies = np.zeros((get_number_of_atoms(res),))
    error_check(lib.tblite_get_result_energies)(
        res, ffi.cast("double*", _energies.ctypes.data)
    )
    return _energies


def get_gradient(res):
    """Retrieve gradient from result container"""
    _gradient = np.zeros((get_number_of_atoms(res), 3))
    error_check(lib.tblite_get_result_gradient)(
        res, ffi.cast("double*", _gradient.ctypes.data)
    )
    return _gradient


def get_virial(res):
    """Retrieve virial from result container"""
    _virial = np.zeros((3, 3))
    error_check(lib.tblite_get_result_virial)(
        res, ffi.cast("double*", _virial.ctypes.data)
    )
    return _virial


def get_charges(res):
    """Retrieve atomic charges from result container"""
    _charges = np.zeros((get_number_of_atoms(res),))
    error_check(lib.tblite_get_result_charges)(
        res, ffi.cast("double*", _charges.ctypes.data)
    )
    return _charges


def get_bond_orders(res):
    """Retrieve Wiberg / Mayer bond orders from result container"""
    _natoms = get_number_of_atoms(res)
    _bond_orders = np.zeros((_natoms, _natoms))
    error_check(lib.tblite_get_result_bond_orders)(
        res, ffi.cast("double*", _bond_orders.ctypes.data)
    )
    return _bond_orders


def get_dipole(res):
    """Retrieve dipole moment from result container"""
    _dipole = np.zeros((3,))
    error_check(lib.tblite_get_result_dipole)(
        res, ffi.cast("double*", _dipole.ctypes.data)
    )
    return _dipole


def get_quadrupole(res):
    """Retrieve quadrupole moment from result container"""
    _quadrupole = np.zeros((6,))
    error_check(lib.tblite_get_result_quadrupole)(
        res, ffi.cast("double*", _quadrupole.ctypes.data)
    )
    return _quadrupole


def get_orbital_energies(res):
    """Retrieve orbital energies from result container"""
    _norb = get_number_of_orbitals(res)
    _emo = np.zeros((_norb,))
    error_check(lib.tblite_get_result_orbital_energies)(
        res, ffi.cast("double*", _emo.ctypes.data)
    )
    return _emo


def get_orbital_occupations(res):
    """Retrieve orbital occupations from result container"""
    _norb = get_number_of_orbitals(res)
    _occ = np.zeros((_norb,))
    error_check(lib.tblite_get_result_orbital_occupations)(
        res, ffi.cast("double*", _occ.ctypes.data)
    )
    return _occ


def _get_ao_matrix(getter):
    """Correctly set allocation for matrix objects before querying the getter"""

    @functools.wraps(getter)
    def with_allocation(res):
        """Get a matrix property from the results object"""
        _norb = get_number_of_orbitals(res)
        _mat = np.zeros((_norb, _norb))
        error_check(getter)(res, ffi.cast("double*", _mat.ctypes.data))
        return _mat

    return with_allocation


def _get_ml_features(getter):
    """Correctly set allocation for matrix objects ml features before querying the getter"""

    @functools.wraps(getter)
    def with_allocation(res):
        """Get a matrix property from the results object"""
        _natoms = ffi.new("int *")
        _nfeatures = ffi.new("int *")
        error_check(lib.tblite_get_result_number_of_atoms)(res, _natoms)
        error_check(lib.tblite_get_result_ml_n_features)(res, _nfeatures)
        _mat = np.zeros((_nfeatures[0], _natoms[0]))
        error_check(getter)(res, ffi.cast("double*", _mat.ctypes.data))
        _mat = _mat.T
        return _mat

    return with_allocation


def _get_w_xtbml(getter):
    """Correctly set allocation for matrix objects w_xtbml before querying the getter"""

    @functools.wraps(getter)
    def with_allocation(res):
        """Get a matrix property from the results object"""
        _natoms = ffi.new("int *")
        error_check(lib.tblite_get_result_number_of_atoms)(res, _natoms)
        _mat = np.zeros((_natoms[0]))
        error_check(getter)(res, ffi.cast("double*", _mat.ctypes.data))
        _mat = _mat
        return _mat

    return with_allocation


def _get_ml_labels(getter):
    """Correctly set allocation for matrix objects ml features before querying the getter"""

    @functools.wraps(getter)
    def with_allocation(res):
        """Get a matrix property from the results object"""
        _nfeatures = ffi.new("int *")
        error_check(lib.tblite_get_result_ml_n_features)(res, _nfeatures)
        labels = list()
        for i in range(1, _nfeatures[0] + 1):
            _index = ffi.new("const int*", i)
            _message = ffi.new("char[]", 512)
            error_check(getter)(res, _message, ffi.NULL, _index)
            labels.append(ffi.string(_message).decode())
        return labels

    return with_allocation


get_orbital_coefficients = _get_ao_matrix(lib.tblite_get_result_orbital_coefficients)
get_density_matrix = _get_ao_matrix(lib.tblite_get_result_density_matrix)
get_overlap_matrix = _get_ao_matrix(lib.tblite_get_result_overlap_matrix)
get_hamiltonian_matrix = _get_ao_matrix(lib.tblite_get_result_hamiltonian_matrix)
get_ml_features = _get_ml_features(lib.tblite_get_result_ml_features)
get_ml_labels = _get_ml_labels(lib.tblite_get_result_ml_labels)


def _delete_calculator(calc) -> None:
    """Delete a tblite calculator object"""
    ptr = ffi.new("tblite_calculator *")
    ptr[0] = calc
    lib.tblite_delete_calculator(ptr)


@context_check
def new_gfn2_calculator(ctx, mol):
    """Create new tblite calculator loaded with GFN2-xTB parametrization data"""
    return ffi.gc(lib.tblite_new_gfn2_calculator(ctx, mol), _delete_calculator)


@context_check
def new_gfn1_calculator(ctx, mol):
    """Create new tblite calculator loaded with GFN1-xTB parametrization data"""
    return ffi.gc(lib.tblite_new_gfn1_calculator(ctx, mol), _delete_calculator)


@context_check
def new_ipea1_calculator(ctx, mol):
    """Create new tblite calculator loaded with IPEA1-xTB parametrization data"""
    return ffi.gc(lib.tblite_new_ipea1_calculator(ctx, mol), _delete_calculator)


@context_check
def new_xtb_calculator(ctx, mol, param):
    """Create new tblite calculator from parametrization records"""
    return ffi.gc(lib.tblite_new_xtb_calculator(ctx, mol, param), _delete_calculator)

def _delete_post_processing(calc) -> None:
    """Delete a tblite pots-processing object"""
    ptr = ffi.new("tblite_post_processing *")
    ptr[0] = calc
    lib.tblite_delete_calculator(ptr)

def new_post_processing(string):
    """Create new tblite calculator from parametrization records"""
    _string = ffi.new("char[]", string.encode("ascii"))
    return ffi.gc(lib.tblite_new_post_processing(_string), _delete_post_processing)

def get_calculator_shell_map(ctx, calc):
    """Retrieve index mapping from shells to atomic centers"""
    _nsh = ffi.new("int *")
    context_check(lib.tblite_get_calculator_shell_count)(ctx, calc, _nsh)
    _map = np.zeros((_nsh[0],), dtype=np.int32)
    context_check(lib.tblite_get_calculator_shell_map)(
        ctx, calc, ffi.cast("int*", _map.ctypes.data)
    )
    return _map


def get_calculator_angular_momenta(ctx, calc):
    """Retrieve angular momenta of shells"""
    _nsh = ffi.new("int *")
    context_check(lib.tblite_get_calculator_shell_count)(ctx, calc, _nsh)
    _am = np.zeros((_nsh[0],), dtype=np.int32)
    context_check(lib.tblite_get_calculator_angular_momenta)(
        ctx, calc, ffi.cast("int*", _am.ctypes.data)
    )
    return _am


def get_calculator_orbital_map(ctx, calc):
    """Retrieve index mapping from atomic orbitals to shells"""
    _nao = ffi.new("int *")
    context_check(lib.tblite_get_calculator_orbital_count)(ctx, calc, _nao)
    _map = np.zeros((_nao[0],), dtype=np.int32)
    context_check(lib.tblite_get_calculator_orbital_map)(
        ctx, calc, ffi.cast("int*", _map.ctypes.data)
    )
    return _map


set_calculator_max_iter = context_check(lib.tblite_set_calculator_max_iter)
set_calculator_accuracy = context_check(lib.tblite_set_calculator_accuracy)
set_calculator_mixer_damping = context_check(lib.tblite_set_calculator_mixer_damping)
set_calculator_guess = context_check(lib.tblite_set_calculator_guess)
set_calculator_temperature = context_check(lib.tblite_set_calculator_temperature)
set_calculator_save_integrals = context_check(lib.tblite_set_calculator_save_integrals)


@context_check
def set_calculator_verbosity(ctx, calc, verbosity: int):
    """Set verbosity in context associated with calculator"""
    lib.tblite_set_context_verbosity(ctx, verbosity)


get_singlepoint = context_check(lib.tblite_get_singlepoint)
get_singlepoint_w_post = context_check(lib.tblite_get_singlepoint_w_post)

def _delete_container(cont) -> None:
    """Delete a tblite container object"""
    ptr = ffi.new("tblite_container *")
    ptr[0] = cont
    lib.tblite_delete_container(ptr)


def new_electric_field(ctx, mol, calc, efield):
    """Create new tblite electric field object"""
    return lib.tblite_new_electric_field(efield)


@context_check
def new_alpb_solvation(ctx, mol, calc, solvent):
    "Create new ALPB solvation model object"
    if isinstance(solvent, str):
        _string = ffi.new("char[]", solvent.encode("ascii"))
        return lib.tblite_new_alpb_solvation_solvent(ctx, mol, calc, _string)
    elif isinstance(solvent, float) or isinstance(solvent, int):
        _eps = float(solvent)
        return lib.tblite_new_alpb_solvation_epsilon(ctx, mol, calc, _eps)
    else:
        raise TBLiteTypeError(
            "Enter desired solvent as string, or enter epsilon value as float or intger."
        )


@context_check
def new_cpcm_solvation(ctx, mol, calc, solvent):
    "Create new ALPB solvation model object"
    if isinstance(solvent, str):
        _string = ffi.new("char[]", solvent.encode("ascii"))
        return lib.tblite_new_cpcm_solvation_solvent(ctx, mol, calc, _string)
    elif isinstance(solvent, float) or isinstance(solvent, int):
        _eps = float(solvent)
        return lib.tblite_new_cpcm_solvation_epsilon(ctx, mol, calc, _eps)
    else:
        raise TBLiteTypeError(
            "Enter desired solvent as string, or enter epsilon value as float or intger."
        )


@context_check
def new_spin_polarization(ctx, mol, calc, wscale: float = 1.0):
    """Create new tblite spin polarization object"""
    return lib.tblite_new_spin_polarization(ctx, mol, calc, wscale)


@context_check
def calculator_push_back(ctx, calc, cont) -> None:
    ptr = ffi.new("tblite_container *")
    ptr[0] = cont
    lib.tblite_calculator_push_back(ctx, calc, ptr)