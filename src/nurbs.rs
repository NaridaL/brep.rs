use crate::V4::V4;
use num_traits::real::Real;
use num_traits::{One, Zero};

#[derive(Debug, Clone)]
pub struct NURBS<T = f64> {
	points: Vec<V4<T>>,
	knots: Vec<T>,
	t_min: T,
	t_max: T,
}

impl<T: Real> NURBS<T> {
	fn new(points: &[V4<T>], degree: u8, knots: &[T], t_min: T, t_max: T) -> NURBS<T> {
		let result = NURBS {
			points: points.to_vec(),
			knots: knots.to_vec(),
			t_min,
			t_max,
		};
		assert_eq!(result.degree(), degree);
		result
	}
	fn degree(&self) -> u8 {
		// knotsLength - points.length - 1 = degree
		(self.knots.len() - self.points.len() - 1) as u8
	}
	fn at4(self, t: T) -> Result<V4<T>, String> {
		let _0: T = Zero::zero();
		let _1: T = One::one();
		let NURBS {
			points,
			knots,
			t_min,
			t_max,
		} = &self;
		let degree = self.degree() as usize;
		assert_between(t, *t_min, *t_max);

		// find s (the spline segment) for the [t] value provided
		let s = self.tInterval(t)?;

		let mut v = points.clone();

		for level in 0..degree {
			// build level l of the pyramid
			for i in (level..degree).rev() {
				let alpha =
					(t - knots[i + s - degree]) / (knots[i + s - level] - knots[i + s - degree]);

				v[i] = v[i - 1] * (_1 - alpha) + v[i] * alpha;
			}
		}

		Ok(v[degree])
	}

	/**
	 * Returns the index of the interval which contains the value t.
	 */
	fn tInterval(&self, t: T) -> Result<usize, String> {
		let knots = &self.knots;
		let degree = self.degree() as usize;
		for s in degree..(knots.len() - 1 - degree) {
			if knots[s] <= t && t <= knots[s + 1] {
				return Ok(s);
			}
		}
		Err("t out of bounds".to_owned())
	}
}

fn assert_between<T: Real>(val: T, min: T, max: T) {
	assert!(min <= val);
	assert!(val <= max);
}
