from numpy.ctypeslib import ndpointer
from ctypes import *
import numpy as np
import os


# Complex number class for the C library
class Complex(Structure):
    _fields_ = [('real', c_double), ('imag', c_double)]


class winterp_interpolator(Structure):
    _fields_ = [('datatype', c_char),
                ('datadim', c_int),
                ('nx', c_int),
                ('ny', c_int),
                ('nz', c_int),
                ('dx', c_double),
                ('dy', c_double),
                ('dz', c_double),
                ('datak', POINTER(Complex))]


# Load the C library
libpath = f"{os.getcwd()+os.sep+os.pardir+os.sep}libwinterp.so"
libWInterp = CDLL(libpath)
libWInterp.connect()

# Define functions from the C library
winterp_getvalue_1d_r = libWInterp.winterp_getvalue_1d_r
winterp_getvalue_1d_r.restype = c_int
winterp_getvalue_1d_r.argtypes = [
    POINTER(winterp_interpolator), c_double, POINTER(c_double)]

winterp_getvalue_1d_c = libWInterp.winterp_getvalue_1d_c
winterp_getvalue_1d_c.restype = c_int
winterp_getvalue_1d_c.argtypes = [
    POINTER(winterp_interpolator), c_double, POINTER(Complex)]

winterp_create_interpolator_1d_r = libWInterp.winterp_create_interpolator_1d_r
winterp_create_interpolator_1d_r.restype = c_int
winterp_create_interpolator_1d_r.argtypes = [
    c_int, c_double, ndpointer(dtype=np.float64), POINTER(winterp_interpolator)]

winterp_create_interpolator_1d_c = libWInterp.winterp_create_interpolator_1d_c
winterp_create_interpolator_1d_c.restype = c_int
winterp_create_interpolator_1d_c.argtypes = [c_int, c_double, ndpointer(
    dtype=[('real', np.float64), ('imag', np.float64)]), POINTER(winterp_interpolator)]


############################
# 2D interpolation functions
############################
winterp_getvalue_2d_r = libWInterp.winterp_getvalue_2d_r
winterp_getvalue_2d_r.restype = c_int
winterp_getvalue_2d_r.argtypes = [
    POINTER(winterp_interpolator), c_double, c_double, POINTER(c_double)]

winterp_getvalue_2d_c = libWInterp.winterp_getvalue_2d_c
winterp_getvalue_2d_c.restype = c_int
winterp_getvalue_2d_c.argtypes = [
    POINTER(winterp_interpolator), c_double, c_double, POINTER(Complex)]

winterp_create_interpolator_2d_r = libWInterp.winterp_create_interpolator_2d_r
winterp_create_interpolator_2d_r.restype = c_int
winterp_create_interpolator_2d_r.argtypes = [c_int, c_int, c_double, c_double, ndpointer(
    dtype=np.float64), POINTER(winterp_interpolator)]

winterp_create_interpolator_2d_c = libWInterp.winterp_create_interpolator_2d_c
winterp_create_interpolator_2d_c.restype = c_int
winterp_create_interpolator_2d_c.argtypes = [c_int, c_int, c_double, c_double, ndpointer(
    dtype=[('real', np.float64), ('imag', np.float64)]), POINTER(winterp_interpolator)]


############################
# 3D interpolation functions
############################
winterp_getvalue_3d_r = libWInterp.winterp_getvalue_3d_r
winterp_getvalue_3d_r.restype = c_int
winterp_getvalue_3d_r.argtypes = [
    POINTER(winterp_interpolator), c_double, c_double, c_double, POINTER(c_double)]

winterp_getvalue_3d_c = libWInterp.winterp_getvalue_3d_c
winterp_getvalue_3d_c.restype = c_int
winterp_getvalue_3d_c.argtypes = [
    POINTER(winterp_interpolator), c_double, c_double, c_double, POINTER(Complex)]

winterp_create_interpolator_3d_r = libWInterp.winterp_create_interpolator_3d_r
winterp_create_interpolator_3d_r.restype = c_int
winterp_create_interpolator_3d_r.argtypes = [c_int, c_int, c_int, c_double, c_double, c_double, ndpointer(
    dtype=np.float64), POINTER(winterp_interpolator)]

winterp_create_interpolator_3d_c = libWInterp.winterp_create_interpolator_3d_c
winterp_create_interpolator_3d_c.restype = c_int
winterp_create_interpolator_3d_c.argtypes = [c_int, c_int, c_int, c_double, c_double, c_double, ndpointer(
    dtype=[('real', np.float64), ('imag', np.float64)]), POINTER(winterp_interpolator)]

# Destroy interp function
winterp_destroy_interpolator = libWInterp.winterp_destroy_interpolator
winterp_destroy_interpolator.restype = c_int
winterp_destroy_interpolator.argtypes = [POINTER(winterp_interpolator)]


###########################
# 1D User Friendly wrappers
###########################
def create_interpolator_1d_r(nx, dx, data):
    interp = winterp_interpolator()
    winterp_create_interpolator_1d_r(nx, dx, data, byref(interp))
    return interp


def create_interpolator_1d_c(nx, dx, data):
    # Create a structured NumPy array with real and imaginary parts
    complex_array = np.zeros(
        nx, dtype=[('real', np.float64), ('imag', np.float64)])
    complex_array['real'] = data.real
    complex_array['imag'] = data.imag

    interp = winterp_interpolator()
    winterp_create_interpolator_1d_c(nx, dx, complex_array, byref(interp))
    return interp


def interpolate_1d_r(interp, x):
    val = c_double()
    winterp_getvalue_1d_r(byref(interp), x, byref(val))
    return val.value


def interpolate_1d_c(interp, x):
    val = Complex()
    winterp_getvalue_1d_c(byref(interp), x, byref(val))
    return val.real + 1j * val.imag


###########################
# 2D User Friendly wrappers
###########################
def create_interpolator_2d_r(nx, ny, dx, dy, data):
    interp = winterp_interpolator()
    winterp_create_interpolator_2d_r(
        nx, ny, dx, dy, data, byref(interp))
    return interp


def create_interpolator_2d_c(nx, ny, dx, dy, data):
    # Create a structured NumPy array with real and imaginary parts
    complex_array = np.zeros((nx, ny), dtype=[(
        'real', np.float64), ('imag', np.float64)])
    complex_array['real'] = data.real
    complex_array['imag'] = data.imag

    interp = winterp_interpolator()
    winterp_create_interpolator_2d_c(
        nx, ny, dx, dy, complex_array, byref(interp))
    return interp


def interpolate_2d_r(interp, x, y):
    val = c_double()
    winterp_getvalue_2d_r(byref(interp), x, y, byref(val))
    return val.value


def interpolate_2d_c(interp, x, y):
    val = Complex()
    winterp_getvalue_2d_c(byref(interp), x, y, byref(val))
    return val.real + 1j * val.imag


###########################
# 3D User Friendly wrappers
###########################
def create_interpolator_3d_r(nx, ny, nz, dx, dy, dz, data):
    interp = winterp_interpolator()
    winterp_create_interpolator_3d_r(
        nx, ny, nz, dx, dy, dz, data, byref(interp))
    return interp


def create_interpolator_3d_c(nx, ny, nz, dx, dy, dz, data):
    # Create a structured NumPy array with real and imaginary parts
    complex_array = np.zeros((nx, ny, nz), dtype=[(
        'real', np.float64), ('imag', np.float64)])
    complex_array['real'] = data.real
    complex_array['imag'] = data.imag

    interp = winterp_interpolator()
    winterp_create_interpolator_3d_c(
        nx, ny, nz, dx, dy, dz, complex_array, byref(interp))
    return interp


def interpolate_3d_r(interp, x, y, z):
    val = c_double()
    winterp_getvalue_3d_r(byref(interp), x, y, z, byref(val))
    return val.value


def interpolate_3d_c(interp, x, y, z):
    val = Complex()
    winterp_getvalue_3d_c(byref(interp), x, y, z, byref(val))
    return val.real + 1j * val.imag


##########################
# Master creation function
##########################
def create_interpolator(data, dx, dy=0.0, dz=0.0):
    data_dim = len(data.shape)
    data_type = ''.join([i for i in str(data.dtype) if not i.isdigit()])
    if data_dim == 1 and data_type == 'float':
        return create_interpolator_1d_r(data.shape[0],
                                        dx,
                                        data)
    elif data_dim == 1 and data_type == 'complex':
        return create_interpolator_1d_c(data.shape[0],
                                        dx,
                                        data)
    elif data_dim == 2 and data_type == 'float':
        return create_interpolator_2d_r(data.shape[0], data.shape[1],
                                        dx, dy,
                                        data)
    elif data_dim == 2 and data_type == 'complex':
        return create_interpolator_2d_c(data.shape[0], data.shape[1],
                                        dx, dy,
                                        data)
    elif data_dim == 3 and data_type == 'float':
        return create_interpolator_3d_r(data.shape[0], data.shape[1], data.shape[2],
                                        dx, dy, dz,
                                        data)
    elif data_dim == 3 and data_type == 'complex':
        return create_interpolator_3d_c(data.shape[0], data.shape[1], data.shape[2],
                                        dx, dy, dz,
                                        data)
    else:
        return None


def interpolate(interpolator, x, y=0, z=0):
    data_dim = interpolator.datadim
    data_type = interpolator.datatype

    if data_dim == 1 and data_type == b'r':
        return interpolate_1d_r(interpolator, x)
    elif data_dim == 1 and data_type == b'c':
        return interpolate_1d_c(interpolator, x)
    elif data_dim == 2 and data_type == b'r':
        return interpolate_2d_r(interpolator, x, y)
    elif data_dim == 2 and data_type == b'c':
        return interpolate_2d_c(interpolator, x, y)
    elif data_dim == 3 and data_type == b'r':
        return interpolate_3d_r(interpolator, x, y, z)
    elif data_dim == 3 and data_type == b'c':
        return interpolate_3d_c(interpolator, x, y, z)
    else:
        return None


def destroy_interpolator(interp):
    winterp_destroy_interpolator(byref(interp))
