use num_traits::real::Real;
use num_traits::{One, Zero};
use std::default::Default;
use std::ops::{Add, Div, Mul, Rem, Sub};

use crate::v4::V4;

#[derive(Debug, Clone)]
struct M4<T: Copy = f64>([T; 16]);
impl<T: Zero + One + Copy> M4<T> {
	fn identity() -> M4<T> {
		let _0: T = Zero::zero();
		let _1: T = One::one();
		M4([
			_1, _0, _0, _0, //
			_0, _1, _0, _0, //
			_0, _0, _1, _0, //
			_0, _0, _0, _1,
		])
	}
}
impl<T: Real> M4<T> {
	fn of_cols(a: V4<T>, b: V4<T>, c: V4<T>, d: V4<T>) -> M4<T> {
		let V4(m00, m10, m20, m30) = a;
		let V4(m01, m11, m21, m31) = b;
		let V4(m02, m12, m22, m32) = c;
		let V4(m03, m13, m23, m33) = d;

		M4([
			m00, m01, m02, m03, //
			m10, m11, m12, m13, //
			m20, m21, m22, m23, //
			m30, m31, m32, m33, //
		])
	}
	fn of_rows(a: V4<T>, b: V4<T>, c: V4<T>, d: V4<T>) -> M4<T> {
		let V4(m00, m01, m02, m03) = a;
		let V4(m10, m11, m12, m13) = b;
		let V4(m20, m21, m22, m23) = c;
		let V4(m30, m31, m32, m33) = d;

		M4([
			m00, m01, m02, m03, //
			m10, m11, m12, m13, //
			m20, m21, m22, m23, //
			m30, m31, m32, m33, //
		])
	}

	fn row(&self, i: usize) -> V4<T> {
		let M4(m) = self;
		V4(m[i * 4], m[i * 4 + 1], m[i * 4 + 2], m[i * 4 + 3])
	}

	fn col(&self, i: usize) -> V4<T> {
		let M4(m) = self;
		V4(m[i], m[i + 4], m[i + 8], m[i + 12])
	}

	fn rows(&self) -> [V4<T>; 4] {
		[self.row(0), self.row(1), self.row(2), self.row(3)]
	}

	/**
	 * Create a rotation matrix for rotating around the Z axis
	 */
	fn rotate_z(radians: T) -> M4<T> {
		let _0: T = Zero::zero();
		let _1: T = One::one();
		let sin = radians.sin();
		let cos = radians.cos();
		M4([
			cos, -sin, _0, _0, //
			sin, cos, _0, _0, //
			_0, _0, _1, _0, //
			_0, _0, _0, _1,
		])
	}
	fn mirror(nx: T, ny: T, nz: T, w: T) -> M4<T> {
		let _0: T = Zero::zero();
		let _1: T = One::one();
		let _2: T = _1 + _1;
		M4([
			_1 - _2 * nx * nx,
			-_2 * ny * nx,
			-_2 * nz * nx,
			_2 * nx * w,
			//
			-_2 * nx * ny,
			_1 - _2 * ny * ny,
			-_2 * nz * ny,
			_2 * ny * w,
			//
			-_2 * nx * nz,
			-_2 * ny * nz,
			_1 - _2 * nz * nz,
			_2 * nz * w,
			//
			_0,
			_0,
			_0,
			_1,
		])
	}

	fn determinant(self) -> T {
		/*
		| a b c d |
		| e f g h |
		| i j k l |
		| m n o p |
		*/
		let M4([a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p]) = self;
		let klop = k * p - l * o;
		let jlnp = j * p - l * n;
		let jkno = j * o - k * n;
		let ilmp = i * p - l * m;
		let ikmo = i * o - k * m;
		let ijmn = i * n - j * m;
		a * (f * klop - g * jlnp + h * jkno) - b * (e * klop - g * ilmp + h * ikmo)
			+ (c * (e * jlnp - f * ilmp + h * ijmn) - d * (e * jkno - f * ikmo + g * ijmn))
	}

	fn inversed(self) -> M4<T> {
		let determinant3 = |a, b, c, d, e, f, g, h, i| -> T {
			a * (e * i - h * h) - b * (d * i - f * g) + c * (d * h - e * g)
		};
		let M4([m00, m01, m02, m03, m10, m11, m12, m13, m20, m21, m22, m23, m30, m31, m32, m33]) =
			self;
		let mut r: [T; 16] = [
			// first compute  cofactor matrix:
			// cofactor of an element is the determinant of the 3x3 matrix gained by removing the column and row belonging
			// to the element
			determinant3(m11, m12, m13, m21, m22, m23, m31, m32, m33),
			determinant3(-m01, m02, m03, m21, m22, m23, m31, m32, m33),
			determinant3(m01, m02, m03, m11, m12, m13, m31, m32, m33),
			determinant3(-m01, m02, m03, m11, m12, m13, m21, m22, m23),
			determinant3(-m10, m12, m13, m20, m22, m23, m30, m32, m33),
			determinant3(m00, m02, m03, m20, m22, m23, m30, m32, m33),
			determinant3(-m00, m02, m03, m10, m12, m13, m30, m32, m33),
			determinant3(m00, m02, m03, m10, m12, m13, m20, m22, m23),
			determinant3(m10, m11, m13, m20, m21, m23, m30, m31, m33),
			determinant3(-m00, m01, m03, m20, m21, m23, m30, m31, m33),
			determinant3(m00, m01, m03, m10, m11, m13, m30, m31, m33),
			determinant3(m00, m01, m03, m10, m11, m13, m20, m21, m23),
			determinant3(-m10, m11, m12, m20, m21, m22, m30, m31, m32),
			determinant3(m00, m01, m02, m20, m21, m22, m30, m31, m32),
			determinant3(-m00, m01, m02, m10, m11, m12, m30, m31, m32),
			determinant3(m00, m01, m02, m10, m11, m12, m20, m21, m22),
		];

		// calculate determinant using laplace expansion (cf https://en.wikipedia.org/wiki/Laplace_expansion),
		// as we already have the cofactors. We multiply a column by a row as the cofactor matrix is transposed.
		let det = m00 * r[0] + m01 * r[4] + m02 * r[8] + m03 * r[12];
		// assert(!isZero(det), 'det may not be zero, i.e. the matrix is not invertible')
		for i in 0..16 {
			r[i] = r[i] / det;
		}
		M4(r)
	}

	fn transposed(self) -> M4<T> {
		let M4([m00, m01, m02, m03, m10, m11, m12, m13, m20, m21, m22, m23, m30, m31, m32, m33]) =
			self;
		M4([
			m00, m10, m20, m30, //
			m01, m11, m21, m31, //
			m02, m12, m22, m32, //
			m03, m13, m23, m33,
		])
	}

	fn frustrum(left: T, right: T, bottom: T, top: T, near: T, far: T) -> M4<T> {
		let _0: T = Zero::zero();
		let _1: T = One::one();
		let _2: T = _1 + _1;

		assert!(_0 < near, "0 < near");
		assert!(near < far, "near < far");

		M4([
			(_2 * near) / (right - left),
			_0,
			(right + left) / (right - left),
			_0,
			//
			_0,
			(_2 * near) / (top - bottom),
			(top + bottom) / (top - bottom),
			_0,
			//
			_0,
			_0,
			-(far + near) / (far - near),
			(-_2 * far * near) / (far - near),
			//
			_0,
			_0,
			-_1,
			_0,
		])
	}
}
macro_rules! multiplex {
    ( $lhs:expr, $rhs:expr, $op:tt ) => { [
		$lhs[0] $op $rhs[0],
		$lhs[1] $op $rhs[1],
		$lhs[2] $op $rhs[2],
		$lhs[3] $op $rhs[3],
		$lhs[4] $op $rhs[4],
		$lhs[5] $op $rhs[5],
		$lhs[6] $op $rhs[6],
		$lhs[7] $op $rhs[7],
		$lhs[8] $op $rhs[8],
		$lhs[9] $op $rhs[9],
		$lhs[10] $op $rhs[10],
		$lhs[11] $op $rhs[11],
		$lhs[12] $op $rhs[12],
		$lhs[13] $op $rhs[13],
		$lhs[14] $op $rhs[14],
		$lhs[15] $op $rhs[15],
	]};
}
macro_rules! multiplex2 {
    ( $lhs:expr, $rhs:expr, $op:tt ) => { [
		$lhs[0] $op $rhs,
		$lhs[1] $op $rhs,
		$lhs[2] $op $rhs,
		$lhs[3] $op $rhs,
		$lhs[4] $op $rhs,
		$lhs[5] $op $rhs,
		$lhs[6] $op $rhs,
		$lhs[7] $op $rhs,
		$lhs[8] $op $rhs,
		$lhs[9] $op $rhs,
		$lhs[10] $op $rhs,
		$lhs[11] $op $rhs,
		$lhs[12] $op $rhs,
		$lhs[13] $op $rhs,
		$lhs[14] $op $rhs,
		$lhs[15] $op $rhs,
	]};
}

impl<RHS: Copy, O: Copy, T: Add<RHS, Output = O> + Copy> Add<M4<RHS>> for M4<T> {
	type Output = M4<O>;
	fn add(self, rhs: M4<RHS>) -> M4<O> {
		M4(multiplex!(self.0, rhs.0, +))
	}
}

impl<RHS: Copy, O: Copy, T: Sub<RHS, Output = O> + Copy> Sub<M4<RHS>> for M4<T> {
	type Output = M4<O>;
	fn sub(self, rhs: M4<RHS>) -> M4<O> {
		M4(multiplex!(self.0, rhs.0, -))
	}
}

impl<RHS: Copy, O: Add<Output = O> + Copy, T: Mul<RHS, Output = O> + Copy> Mul<M4<RHS>> for M4<T> {
	type Output = M4<O>;
	fn mul(self, rhs: M4<RHS>) -> M4<O> {
		let a = self.0;
		let b = rhs.0;
		M4([
			a[0] * b[0] + a[1] * b[4] + (a[2] * b[8] + a[3] * b[12]),
			a[0] * b[1] + a[1] * b[5] + (a[2] * b[9] + a[3] * b[13]),
			a[0] * b[2] + a[1] * b[6] + (a[2] * b[10] + a[3] * b[14]),
			a[0] * b[3] + a[1] * b[7] + (a[2] * b[11] + a[3] * b[15]),
			//
			a[4] * b[0] + a[5] * b[4] + (a[6] * b[8] + a[7] * b[12]),
			a[4] * b[1] + a[5] * b[5] + (a[6] * b[9] + a[7] * b[13]),
			a[4] * b[2] + a[5] * b[6] + (a[6] * b[10] + a[7] * b[14]),
			a[4] * b[3] + a[5] * b[7] + (a[6] * b[11] + a[7] * b[15]),
			//
			a[8] * b[0] + a[9] * b[4] + (a[10] * b[8] + a[11] * b[12]),
			a[8] * b[1] + a[9] * b[5] + (a[10] * b[9] + a[11] * b[13]),
			a[8] * b[2] + a[9] * b[6] + (a[10] * b[10] + a[11] * b[14]),
			a[8] * b[3] + a[9] * b[7] + (a[10] * b[11] + a[11] * b[15]),
			//
			a[12] * b[0] + a[13] * b[4] + (a[14] * b[8] + a[15] * b[12]),
			a[12] * b[1] + a[13] * b[5] + (a[14] * b[9] + a[15] * b[13]),
			a[12] * b[2] + a[13] * b[6] + (a[14] * b[10] + a[15] * b[14]),
			a[12] * b[3] + a[13] * b[7] + (a[14] * b[11] + a[15] * b[15]),
		])
	}
}
impl<RHS: Copy, O: Add<Output = O> + Copy, T: Mul<RHS, Output = O> + Copy> Mul<V4<RHS>> for M4<T> {
	type Output = V4<O>;
	fn mul(self, rhs: V4<RHS>) -> V4<O> {
		let M4([m00, m01, m02, m03, m10, m11, m12, m13, m20, m21, m22, m23, m30, m31, m32, m33]) =
			self;
		let V4(x, y, z, w) = rhs;

		V4(
			m00 * x + m01 * y + m02 * z + m03 * w,
			m10 * x + m11 * y + m12 * z + m13 * w,
			m20 * x + m21 * y + m22 * z + m23 * w,
			m30 * x + m31 * y + m32 * z + m33 * w,
		)
	}
}

macro_rules! multiplex2_impl {
    ( $tt:ty ) => {
		impl<O: Copy, T: Mul<$tt, Output = O> + Copy> Mul<$tt> for M4<T> {
			type Output = M4<O>;
			fn mul(self, rhs: $tt) -> M4<O> {
				M4(multiplex2!(self.0, rhs, *))
			}
		}
		impl<O: Copy, T: Div<$tt, Output = O> + Copy> Div<$tt> for M4<T> {
			type Output = M4<O>;
			fn div(self, rhs: $tt) -> M4<O> {
				M4(multiplex2!(self.0, rhs, /))
			}
		}
	}
}
multiplex2_impl!(f32);
multiplex2_impl!(f64);
multiplex2_impl!(i32);
multiplex2_impl!(i64);
