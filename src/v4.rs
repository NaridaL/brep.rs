use num_traits::real::Real;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Mul;
use std::ops::Sub;

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct V4<T = f64>(pub T, pub T, pub T, pub T);
// {
//     x: T,
//     y: T,
//     z: T,
//     w: T,
// }
impl<T> V4<T> {
	pub fn new(x: T, y: T, z: T, w: T) -> V4<T> {
		V4(x, y, z, w)
	}
}

impl<T: Add<Output = T>> Add<V4<T>> for V4<T> {
	type Output = V4<T>;
	fn add(self, rhs: V4<T>) -> V4<T> {
		let V4(x, y, z, w) = self;
		let V4(bx, by, bz, bw) = rhs;
		V4(x + bx, y + by, z + bz, w + bw)
	}
}
impl<T: Sub<Output = T>> Sub<V4<T>> for V4<T> {
	type Output = V4<T>;
	fn sub(self, rhs: V4<T>) -> V4<T> {
		let V4(x, y, z, w) = self;
		let V4(bx, by, bz, bw) = rhs;
		V4(x - bx, y - by, z - bz, w - bw)
	}
}

impl<T: AddAssign> AddAssign<V4<T>> for V4<T> {
	fn add_assign(&mut self, _rhs: V4<T>) {
		let V4(ref mut x1, ref mut y1, ref mut z1, ref mut w1) = self;
		let V4(x2, y2, z2, w2) = _rhs;
		*x1 += x2;
		*y1 += y2;
		*z1 += z2;
		*w1 += w2;
	}
}

impl<T: Real> Mul<T> for V4<T> {
	type Output = T;
	fn mul(self, _rhs: T) -> T {
		let V4(x1, y1, z1, w1) = self;

		x1 * _rhs + y1 * _rhs + (z1 * _rhs + w1 * _rhs)
	}
}
impl<T: Mul<Output = T> + std::marker::Copy> Mul<V4<T>> for V4<T> {
	type Output = V4<T>;
	fn mul(self, _rhs: V4<T>) -> V4<T> {
		let V4(x1, y1, z1, w1) = self;
		let V4(x2, y2, z2, w2) = _rhs;

		V4(x1 * x2, y1 * y2, z1 * z2, w1 * w2)
	}
}
