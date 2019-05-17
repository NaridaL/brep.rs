use bit_vec::BitVec;
use itertools::Itertools;
use num_traits::real::Real;
use ord_subset::OrdSubsetIterExt;
use std::fs::File;
use std::io::Write;
use std::iter::FromIterator;

pub fn earcut(data: &[P2], hole_indices: &[usize], dim: u8) {
	let has_holes = !hole_indices.is_empty();
	let e = |i: usize| data[i];
	let nv = data.len();
	let mut data2 = Vec::with_capacity(nv + hole_indices.len() + 2);
	data2.extend_from_slice(data);
	data2[0..nv].copy_from_slice(data);

	let mut next_unused_index = Vec::with_capacity(data2.len());
	for i in 0..nv {
		next_unused_index.push(i + 1);
	}
	let outer_loop_end = hole_indices.get(0).cloned().unwrap_or(data.len());
	next_unused_index[outer_loop_end - 1] = 0;
	for (i, &start) in hole_indices.iter().enumerate() {
		let end = hole_indices.get(i + 1).cloned().unwrap_or(data.len());
		next_unused_index[end - 1] = start;
	}

	struct TriangleCycleIterator<'a, T: Copy> {
		w: (T, T, T),
		i: std::slice::Iter<'a, T>,
	}
	impl<'a, T: Copy> TriangleCycleIterator<'a, T> {
		fn new(data: &[T]) -> TriangleCycleIterator<T> {
			TriangleCycleIterator {
				w: (
					data[data.len() - 3],
					data[data.len() - 2],
					data[data.len() - 1],
				),
				i: data.iter(),
			}
		}
	}
	impl<'a, T: Copy> Iterator for TriangleCycleIterator<'a, T> {
		type Item = (T, T, T);
		fn next(&mut self) -> Option<(T, T, T)> {
			self.i.next().map(|next| {
				self.w.0 = self.w.1;
				self.w.1 = self.w.2;
				self.w.2 = *next;

				self.w
			})
		}
	}

	let mut is_convex = BitVec::with_capacity(data2.len());
	for i in 0..outer_loop_end {
		is_convex.push(triangle_is_convex(
			e((i + outer_loop_end - 1) % outer_loop_end),
			e(i),
			e((i + 1) % outer_loop_end),
		))
	}

	for (i, &start) in hole_indices.iter().enumerate() {
		let end = hole_indices.get(i + 1).cloned().unwrap_or(data.len());
		next_unused_index[end - 1] = start;
		let hole_data = &data[start..end];
		let hole_i = start
			+ hole_data
				.iter()
				.enumerate()
				.ord_subset_min_by_key(|(_, (x, _))| x)
				.unwrap()
				.0;
		let hole_i_prev = start + ((hole_i - start) + 1) % (end - start);
		let (bridge_i_prev, bridge_i) =
			find_hole_bridge(&next_unused_index, &is_convex, hole_i, e).unwrap();

		// bridge_i_prev -> bridge_i -> bridge_i_next and
		// hole_i_prev -> hole_i -> hole_i_next
		// become
		// bridge_i_prev -> bridge_i -> hole_i2 -> hole_i_next
		// hole_i_prev -> hole_i -> bridge_i2 -> bridge_i_next
		let bridge_i_next = next_unused_index[bridge_i];
		let hole_i_next = next_unused_index[hole_i];

		let hole_i2 = data2.len();
		next_unused_index[bridge_i] = hole_i2;
		data2.push(e(hole_i));
		next_unused_index.push(hole_i_next);

		let bridge_i2 = data2.len();
		next_unused_index[hole_i] = bridge_i2;
		data2.push(e(bridge_i));
		next_unused_index.push(bridge_i_next);

		is_convex.set(
			bridge_i,
			triangle_is_convex(e(bridge_i_prev), e(bridge_i), e(hole_i2)),
		);
		is_convex.set(
			hole_i,
			triangle_is_convex(e(hole_i_prev), e(hole_i), e(bridge_i2)),
		);
		is_convex.push(triangle_is_convex(e(bridge_i), e(hole_i2), e(hole_i_next)));
		is_convex.push(triangle_is_convex(
			e(hole_i),
			e(bridge_i2),
			e(bridge_i_next),
		));
	}
}

fn archimedean(steps: usize, start: f64, turns: f64) -> impl Iterator<Item = (f64, f64)> {
	assert_ne!(steps, 0, "Steps should be at least 1.");
	let turns_per_step = (turns - start) / steps as f64;
	(0..=steps).map(move |i| {
		let rad = start + (i as f64 * turns_per_step);
		(rad * rad.cos(), rad * rad.abs().sin())
	})
}

static PI: f64 = std::f64::consts::PI;
static TAU: f64 = PI * 2.0;

fn apply(p: P2, m: (f64, f64, f64, f64)) -> P2 {
	let (x, y) = p;
	let (a, b, c, d) = m;
	(a * x + b * y, c * x + d * y)
}

fn multi_archi(center_distance: f64, num: usize) -> impl Iterator<Item = (f64, f64)> {
	// const center_distance: f64 = 14.0;
	// const num: usize = 5;
	fn rot(x: f64) -> (f64, f64, f64, f64) {
		(x.cos(), -x.sin(), x.sin(), x.cos())
	}

	let closures = (0..num).map(move |i| {
		move |(x, y)| {
			let trans = (x + center_distance, y);
			let mat = rot(TAU * i as f64 / num as f64);
			apply(trans, mat)
		}
	});
	closures.flat_map(|c| archimedean(128, -7.0, 9.0).map(c))
}

fn earcut_no_holes(data: &[f64]) -> Vec<usize> {
	let modulo = |a: usize, b: usize| -> usize { ((a % b) + b) % b };
	let nv = data.len() / 2;
	let e = |i: usize| (data[i * 2], data[i * 2 + 1]);
	let mut triangles: Vec<usize> = Vec::with_capacity(3 * (nv - 2));
	let mut is_convex = BitVec::from_fn(nv, |i| {
		triangle_is_convex(e((i + nv - 1) % nv), e(i), e((i + 1) % nv))
	});
	// println!("is_convex {:?}", is_convex);
	let mut is_used = BitVec::from_elem(nv, false);
	let mut vertices_left = nv;
	let mut current_indexes = [0, 1, 2, 3 % nv];
	let mut next_unused_index = Vec::from_iter((0..nv).map(|i| (i + 1) % nv));
	fn unused_index_after(next_unused_index: &[usize], i: usize) -> usize {
		next_unused_index[i]
	};
	// fn unused_index_after(is_used: &BitVec, mut i: usize) -> usize {
	// 	next_unused_index
	// 	loop {
	// 		i = (i + 1) % is_used.len();
	// 		if !is_used[i] {
	// 			break;
	// 		}
	// 	}
	// 	i
	// };

	let mut min_x = f64::max_value();
	let mut max_x = f64::min_value();
	let mut min_y = f64::max_value();
	let mut max_y = f64::min_value();
	for (&x, &y) in data.iter().tuples() {
		min_x = min_x.min(x);
		max_x = max_x.max(x);
		min_y = min_y.min(y);
		max_y = max_y.max(y);
	}
	let width_inv = 1.0 / (max_x - min_x);
	let height_inv = 1.0 / (max_y - min_y);

	fn lerp_inv(a: f64, b: f64, t: f64) -> f64 {
		(t - a) / (b - a)
	}

	let z_order_val = |(x, y)| {
		interleave(
			(u16::max_value() as f64 * (x - min_x) * width_inv) as u16,
			(u16::max_value() as f64 * (y - min_y) * height_inv) as u16,
		)
	};
	let mut z_order_values =
		Vec::from_iter(data.iter().cloned().tuples().map(z_order_val).enumerate());
	z_order_values.sort_by_key(|(_, z): &(_, u32)| *z);
	let z_index_arr = {
		let mut x = vec![0; nv];
		for (i, j) in z_order_values.iter().map(|(j, _)| *j).enumerate() {
			x[j] = i;
		}
		x
	};

	fn is_ear<F>(
		is_convex: &BitVec,
		is_used: &BitVec,
		e: impl Fn(usize) -> (f64, f64),
		i: usize,
		j: usize,
		k: usize,
	) -> bool {
		// println!("is_convex {:?}", is_convex[j]);
		let (ax, ay) = e(i);
		let (bx, by) = e(j);
		let (cx, cy) = e(k);
		is_convex[j]
			&& is_used.iter().enumerate().all(|(vi, used)| {
				used || vi == i
					|| vi == j || vi == k
					|| !triangle_contains_point(ax, ay, bx, by, cx, cy, e(vi))
			})
	};

	fn is_ear_zorder(
		is_convex: &BitVec,
		is_used: &BitVec,
		z_order_values: &[(usize, u32)],
		z_order_val: impl Fn((f64, f64)) -> u32,
		z_index_arr: &[usize],
		e: impl Fn(usize) -> (f64, f64),
		i: usize,
		j: usize,
		k: usize,
	) -> bool {
		if !is_convex[j] {
			return false;
		}

		// println!("is_ear_zorder {} {} {}", i, j, k);

		let nv = z_order_values.len();

		let (ax, ay) = e(i);
		let (bx, by) = e(j);
		let (cx, cy) = e(k);

		let min_triang_x = ax.min(bx).min(cx);
		let max_triang_x = ax.max(bx).max(cx);
		let min_triang_y = ay.min(by).min(cy);
		let max_triang_y = ay.max(by).max(cy);

		let min_triang_z = z_order_val((min_triang_x, min_triang_y));
		let max_triang_z = z_order_val((max_triang_x, max_triang_y));
		let min_triang_z_x = min_triang_z & 0x5555_5555;
		let min_triang_z_y = min_triang_z & 0xAAAA_AAAA;
		let max_triang_z_x = max_triang_z & 0x5555_5555;
		let max_triang_z_y = max_triang_z & 0xAAAA_AAAA;
		// let min_triang_z = 0;
		// let max_triang_z = u32::max_value();

		let a = || {
			(0..z_index_arr[j])
				.rev()
				.map(|zi| z_order_values[zi])
				.take_while(|&(_, z)| min_triang_z <= z)
		};

		let b = || {
			((z_index_arr[j] + 1)..nv)
				.map(|zi| z_order_values[zi])
				.take_while(|&(_, z)| z <= max_triang_z)
		};

		// let correct_answer = (0..nv).all(|vi| {
		// 	is_used[vi]
		// 		|| vi == i || vi == j
		// 		|| vi == k || !triangle_contains_point(ax, ay, bx, by, cx, cy, e(vi * 2), e(vi * 2 + 1))
		// });
		// let mut tested = Vec::new();
		// (0..nv).all(|vi| {
		let result = a().interleave(b()).all(|(vi, z)| {
			let z_order_x_part = z & 0x5555_5555;
			let z_order_y_part = z & 0xAAAA_AAAA;
			// tested.push(vi);
			is_used[vi]
				|| vi == i || vi == j
				|| vi == k || z_order_x_part > max_triang_z_x
				|| z_order_x_part < min_triang_z_x
				|| z_order_y_part > max_triang_z_y
				|| z_order_y_part < min_triang_z_y
				|| !triangle_contains_point(ax, ay, bx, by, cx, cy, e(vi))
		});
		// if correct_answer != result {
		// 	println!(
		// 		"triangle bounds: {} {} {} {} z  {} {}",
		// 		min_triang_x, min_triang_y, max_triang_x, max_triang_y, min_triang_z, max_triang_z
		// 	);
		// 	a().interleave(b())
		// 		.for_each(|(vi, z)| println!("vi, z {} {}", vi, z));
		// 	let fail_vi = (0..nv)
		// 		.find(|vi| {
		// 			!tested.contains(vi)
		// 				&& triangle_contains_point(ax, ay, bx, by, cx, cy, e(vi * 2), e(vi * 2 + 1))
		// 		})
		// 		.unwrap();
		// 	let (_, fail_vi_z) = z_order_values[z_index_arr[fail_vi]];
		// 	println!(
		// 		"fail_vi_z {} {}     {} {}",
		// 		fail_vi,
		// 		fail_vi_z,
		// 		e(fail_vi * 2),
		// 		e(fail_vi * 2 + 1)
		// 	);
		// }
		// assert_eq!(correct_answer, result);
		result

		// look for points inside the triangle in both directions
		// loop {
		// 	let (_,p_z) = z_order_values[p_zi];
		// 	_
		// 	while (p && p.z >= minZ && n && n.z <= maxZ)
		//     if (p !== ear.prev && p !== ear.next &&
		//         pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) &&
		//         area(p.prev, p, p.next) >= 0) return false;
		//     p = p.prevZ;

		//     if (n !== ear.prev && n !== ear.next &&
		//         pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, n.x, n.y) &&
		//         area(n.prev, n, n.next) >= 0) return false;
		//     n = n.nextZ;
		// }

		// // look for remaining points in decreasing z-order
		// while (p && p.z >= minZ) {
		//     if (p !== ear.prev && p !== ear.next &&
		//         pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) &&
		//         area(p.prev, p, p.next) >= 0) return false;
		//     p = p.prevZ;
		// }

		// // look for remaining points in increasing z-order
		// while (n && n.z <= maxZ) {
		//     if (n !== ear.prev && n !== ear.next &&
		//         pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, n.x, n.y) &&
		//         area(n.prev, n, n.next) >= 0) return false;
		//     n = n.nextZ;
		// }

		// return true;

		// println!("is_convex {:?}", is_convex[j]);
		// is_used.iter().enumerate().all(|(vi, used)| {
		// 	used || vi == i
		// 		|| vi == j || vi == k
		// 		|| !triangle_contains_point(
		// 			e(i * 2),
		// 			e(i * 2 + 1),
		// 			e(j * 2),
		// 			e(j * 2 + 1),
		// 			e(k * 2),
		// 			e(k * 2 + 1),
		// 			e(vi * 2),
		// 			e(vi * 2 + 1),
		// 		)
		// })
	};

	while vertices_left >= 3 {
		{
			let [_, i, j, k] = current_indexes;
			// println!("{:?}", current_indexes);
			let start_end_index = k;
			// find next ear
			while !is_ear_zorder(
				// while !is_ear(
				&is_convex,
				&is_used,
				&z_order_values,
				z_order_val,
				&z_index_arr,
				e,
				current_indexes[1],
				current_indexes[2],
				current_indexes[3],
			) {
				current_indexes[0] = current_indexes[1];
				current_indexes[1] = current_indexes[2];
				current_indexes[2] = current_indexes[3];
				current_indexes[3] = unused_index_after(&next_unused_index, current_indexes[2]);
				// shift_curent_indexes()
				// println!("{:?}", current_indexes);
				if start_end_index == current_indexes[3] {
					println!("endless loop");
					is_used.iter().enumerate().for_each(|(vi, used)| {
						if !used {
							println!("{:?}", e(vi));
						}
					});
					return triangles;
				}
			}
		}
		{
			let [_, i, j, k] = current_indexes;
			triangles.push(i);
			triangles.push(j);
			triangles.push(k);
			is_used.set(j, true);
			vertices_left = vertices_left - 1;
			// println!("added, left: {}", vertices_left);
		}
		// remove j from the current window
		next_unused_index[current_indexes[1]] = next_unused_index[current_indexes[2]];
		current_indexes[2] = current_indexes[3];
		current_indexes[3] = unused_index_after(&next_unused_index, current_indexes[2]);
		{
			let [h, i, j, k] = current_indexes;
			// fix convex values at i and j (previously k)
			is_convex.set(i, is_convex[i] || triangle_is_convex(e(h), e(i), e(j)));
			is_convex.set(j, is_convex[j] || triangle_is_convex(e(i), e(j), e(k)));
		}
		// current_indexes[0] = current_indexes[1];
		// current_indexes[1] = current_indexes[2];
		// current_indexes[2] = current_indexes[3];
		// current_indexes[3] = unused_index_after(&next_unused_index, current_indexes[2]);
	}
	triangles
}

fn double_area((ax, ay): P2, (bx, by): P2, (cx, cy): P2) -> f64 {
	(bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
}

fn triangle_is_convex(a: P2, b: P2, c: P2) -> bool {
	// println!("{} {} {} {} {} {} ", ax, ay, bx, by, cx, cy);
	// println!("{}", (bx - ax) * (cy - ay) - (by - ay) * (cx - ax));
	double_area(a, b, c) > 0.0
}

fn triangle_contains_point(
	ax: f64,
	ay: f64,
	bx: f64,
	by: f64,
	cx: f64,
	cy: f64,
	(px, py): P2,
) -> bool {
	// println!(
	// 	"triangle_contains_point {} {} {} {} {} {}  {} {} ",
	// 	ax, ay, bx, by, cx, cy, px, py
	// );
	let result = (cx - px) * (ay - py) - (ax - px) * (cy - py) >= 0.0
		&& (ax - px) * (by - py) - (bx - px) * (ay - py) >= 0.0
		&& (bx - px) * (cy - py) - (cx - px) * (by - py) >= 0.0;
	// println!("{}", result);
	result
}
fn foo(x: f64, y: f64, a: f64, c: f64) -> f64 {
	x.powi(2) + y.powi(2)
		- a * ((x - y).sin().powi(3).exp() + (-x - y).sin().powi(3).exp()).powi(2)
		- c
}

fn ddd(f: impl Fn(f64, f64) -> f64, v: P2) -> (f64, f64) {
	let eps = 1e-6;
	let (x, y) = v;
	let fxy = f(x, y);
	((f(x + eps, y) - fxy) / eps, (f(x, y + eps) - fxy) / eps)
}

fn curve_point(f: &impl Fn(f64, f64) -> f64, start: (f64, f64)) -> (f64, f64) {
	let mut p = start;
	for _ in 0..8 {
		let (px, py) = p;
		let fp = f(px, py);
		let (dfpdx, dfpdy) = ddd(f, p);
		let scale = fp / (dfpdx * dfpdx + dfpdy * dfpdy);
		p = (px - scale * dfpdx, py - scale * dfpdy);
		// println!("p {:?}       f {:?}", p, f(px, py));
	}
	p
}

fn distance(a: P2, b: P2) -> f64 {
	let (ax, ay) = a;
	let (bx, by) = b;
	length((ax - bx, ay - by))
}

type P2 = (f64, f64);
fn to_length(v: P2, nl: f64) -> P2 {
	let l = length(v);
	let (x, y) = v;
	(x * nl / l, y * nl / l)
}
fn length(v: P2) -> f64 {
	let (x, y) = v;
	(x * x + y * y).sqrt()
}

fn bit_space_out(a: u16) -> u32 {
	let mut x: u32 = a.into();
	x = (x | (x << 8)) & 0x00FF00FF;
	x = (x | (x << 4)) & 0x0F0F0F0F;
	x = (x | (x << 2)) & 0x33333333;
	x = (x | (x << 1)) & 0x55555555;
	x
}

struct AABB2 {
	min: P2,
	max: P2,
}
static INF: f64 = std::f64::INFINITY;
impl AABB2 {
	fn new() -> AABB2 {
		AABB2 {
			min: (INF, INF),
			max: (-INF, -INF),
		}
	}
	fn add_iter(self: &mut Self, iter: impl IntoIterator<Item = P2>) {
		for (x, y) in iter {
			self.min.0 = self.min.0.min(x);
			self.max.0 = self.max.0.max(x);

			self.min.1 = self.min.1.min(y);
			self.max.1 = self.max.1.max(y);
		}
	}
}

fn interleave(a: u16, b: u16) -> u32 {
	bit_space_out(a) | (bit_space_out(b) << 1)
}

// fn cycle_windows<T: Copy, K: PartialEq>(
// 	into_iter: impl IntoIterator<Item = T>,
// 	key: impl Fn(T) -> K,
// ) -> impl Iterator<Item = (T, T, T)> {
// 	let mut iter = into_iter.into_iter();
// 	let stop: T = iter.next().unwrap();
// 	let stop_key = key(stop);
// 	let mut a = stop;
// 	let mut b = iter.next().unwrap();
// 	let mut c = iter.next().unwrap();
// 	std::iter::from_fn(move || {
// 		let tuple = (a, b, c);
// 			a = b;
// 			b = c;
// 			c = iter.next().unwrap();

// 		if key(b) == stop_key {
// 			None
// 		} else {
// 			Some((a, b, c))
// 		}
// 	})
// }
fn cycle_windows<T: Copy, K: PartialEq, R>(
	into_iter: impl IntoIterator<Item = T>,
	key: impl Fn(T) -> K,
	mut f: impl FnMut([T; 3]) -> Option<R>,
) -> Option<R> {
	let mut iter = into_iter.into_iter();

	let mut a = iter.next().unwrap();
	let mut b = iter.next().unwrap();
	let mut c = iter.next().unwrap();

	let stop = a;
	let stop_key = key(stop);
	loop {
		if let Some(x) = f([a, b, c]) {
			return Some(x);
		};

		a = b;
		b = c;
		c = iter.next().unwrap();
		if key(a) == stop_key {
			return None;
		}
	}
}

// David Eberly's algorithm for finding a bridge between hole and outer polygon
// return (i_prev, i), where i the the bridge point, and i_prev is the point preceding
// the bridge point
fn find_hole_bridge(
	next_unused_index: &[usize],
	is_convex: &BitVec,
	hole_left_most: usize,
	e: impl Fn(usize) -> (f64, f64),
) -> Option<(usize, usize)> {
	let (hole_point_x, hole_point_y) = e(hole_left_most); // hole x/y
	let mut pi = 0;
	// x at which the ray from hole_left_most in -x direction intersects nearest segment
	let mut qx = -std::f64::INFINITY;
	let mut m = None;

	// find a segment intersected by a ray from the hole's leftmost point to the left;
	// segment's endpoint with lesser x will be potential connection point
	if let Some(xxx) = cycle_windows(
		std::iter::repeat_with(|| {
			pi = next_unused_index[pi];
			pi
		})
		.map(|i| (i, e(i))),
		|(i, _)| i,
		|[(i_prev, _), (i, (px, py)), (i_next, (p_next_x, p_next_y))]| {
			if p_next_y <= hole_point_y && hole_point_y <= py && p_next_y != py {
				let x = px + (hole_point_y - py) * (p_next_x - px) / (p_next_y - py);
				if x <= hole_point_x && x > qx {
					qx = x;
					if x == hole_point_x {
						return Some(if hole_point_y == py {
							(i_prev, i)
						} else {
							(i, i_next)
						});
					}
					m = Some(if px < p_next_x {
						(i_prev, i)
					} else {
						(i, i_next)
					})
				}
			}
			None
		},
	) {
		return Some(xxx);
	};
	m.map(|mut m| {
		// look for points inside the triangle of hole point, segment intersection and endpoint;
		// if there are no points found, we have a valid connection;
		// otherwise choose the point of the minimum angle with the ray as connection point

		let stop = m;
		let mut tan_min = INF;

		// p = m.next;
		let (mut mi_prev, mut mi) = m;
		let (mut pi_prev, mut pi) = m;
		let (mx, my) = e(mi);

		loop {
			pi_prev = pi;
			pi = next_unused_index[pi];

			if pi == mi {
				break;
			}

			let (px, py) = e(pi);
			if hole_point_x >= px
				&& px >= mx && hole_point_x != px
				&& triangle_contains_point(
					if hole_point_y < my { hole_point_x } else { qx },
					hole_point_y,
					mx,
					my,
					if hole_point_y < my { qx } else { hole_point_x },
					hole_point_y,
					(px, py),
				) {
				let tan = (hole_point_y - py).abs() / (hole_point_x - px); // tangential

				if (tan < tan_min || (tan == tan_min && px > mx)) && !is_convex[pi] {
					mi = pi;
					mi_prev = pi_prev;
					tan_min = tan;
				}
			}
		}

		(mi_prev, mi)
	})
}
type M2 = (f64, f64, f64, f64);
fn mm((a, b, c, d): M2, (e, f, g, h): M2) -> M2 {
	(a * e + b * g, a * f + b * h, c * e + d * g, c * f + d * h)
}

fn menger_wipe(iterations: u32) -> (Vec<P2>, Vec<usize>) {
	fn fn_menger_recursive(
		points: &mut Vec<P2>,
		hole_starts: &mut Vec<usize>,
		lvl: u32,
		(ox, oy): P2,
		s: f64,
	) {
		println!("{:?} {:?}", (ox, oy), s);
		hole_starts.push(points.len());
		for (x, y) in &[
			(1.0 / 3.0, 1.0 / 3.0),
			(1.0 / 3.0, 2.0 / 3.0),
			(2.0 / 3.0, 2.0 / 3.0),
			(2.0 / 3.0, 1.0 / 3.0),
		] {
			points.push((ox + s * x, oy + s * y));
		}

		if 0 != lvl {
			for i in 0..3 {
				for j in 0..3 {
					if j != 1 || i != 1 {
						println!("i {} j {}", i, j);
						fn_menger_recursive(
							points,
							hole_starts,
							lvl - 1,
							(ox + (i as f64 / 3.0) * s, oy + (j as f64 / 3.0) * s),
							s / 3.0,
						);
					}
				}
			}
		}
	};

	let hole_count = (8usize.pow(1 + iterations) - 1) / 7;
	let mut points = Vec::with_capacity((1 + hole_count) * 4);
	let mut hole_starts = Vec::with_capacity(hole_count);
	points.push((0.0, 0.0));
	points.push((1.0, 0.0));
	points.push((1.0, 1.0));
	points.push((0.0, 1.0));
	fn_menger_recursive(&mut points, &mut hole_starts, iterations, (0.0, 0.0), 1.0);
	(points, hole_starts)
}

// check if a polygon diagonal is locally inside the polygon
fn locallyInside(a: P2, b: P2, c: P2, k: P2) -> bool {
	if double_area(a, b, c) < 0.0 {
		double_area(b, k, c) >= 0.0 && double_area(b, a, k) >= 0.0
	} else {
		double_area(b, k, a) < 0.0 || double_area(b, c, k) < 0.0
	}
}

fn followAlgorithm2d(
	ic: &impl Fn(f64, f64) -> f64,
	start_p: (f64, f64),
	step_length: f64,
) -> (Vec<(f64, f64)>, Vec<(f64, f64)>) {
	let (ddx, ddy) = ddd(ic, start_p);
	let start_tangent = to_length((-ddy, ddx), step_length);
	let mut i = 0;
	let mut p = start_p;
	let mut tangent = start_tangent;
	let fullLoop = false;
	let mut points: Vec<P2> = Vec::new();
	let mut tangents: Vec<P2> = Vec::new();
	loop {
		// println!("tangent {:?}", start_tangent);
		points.push(p);
		tangents.push(tangent);
		let (px, py) = p;
		let (tx, ty) = tangent;
		let search_start = (px + tx, py + ty);
		let new_p = curve_point(&ic, search_start);
		let (ddx, ddy) = ddd(ic, new_p);
		let new_tangent = to_length((-ddy, ddx), step_length);
		//const reversedDir = p.minus(prevp).dot(tangent) < 0
		// assert(!p.equals(newP))
		// check if we passed a singularity
		// if (tangent.dot(newTangent) < 0) {
		// 	const singularity = newtonIterate2d(ic.x, ic.y, p.x, p.y)!
		// 	if (eq0(ic(singularity.x, singularity.y)) && singularity.distanceTo(p) < abs(stepLength)) {
		// 		// end on this point
		// 		points.push(singularity)
		// 		tangents.push(p.to(singularity))
		// 		break
		// 	} else {
		// 		throw new Error()
		// 	}
		// }
		// check for endP
		// if (endP && p.equals(endP)) {
		// 	break
		// }
		// check if loop
		if i > 4 && distance(p, start_p) <= step_length.abs() {
			break;
		}
		// check if out of bounds
		// assert(eq0(ic(newP.x, newP.y), NLA_PRECISION * 2), p, newP, searchStart, ic(newP.x, newP.y))
		tangent = new_tangent;
		p = new_p;
		i = i + 1;
		if i == 100000 {
			break;
		};
	}
	// assert!(i < 1000);

	//assert(points.length > 6)
	(points, tangents)
}

use std::iter::once;

fn pair_iter<T>(pair: (T, T)) -> impl Iterator<Item = T> {
	let (a, b) = pair;
	once(a).chain(once(b))
}

#[macro_use]
#[cfg(test)]
mod test {
	extern crate time_test;

	use super::earcut_no_holes;
	use super::foo;
	use super::lerp;
	use super::write_svg;
	use num_traits::float::FloatConst;
	use std::iter::once;

	#[test]
	fn _3vertices() {
		assert_eq!(earcut_no_holes(&[0.0, 0.0, 1.0, 0.0, 0.0, 1.0,]), [1, 2, 0]);
	}

	#[test]
	fn _4vertices() {
		assert_eq!(
			earcut_no_holes(&[0.0, 0.0, 1.0, 0.0, 0.5, 0.5, 0.0, 1.0,]),
			[2, 3, 0, 2, 0, 1]
		);
		assert_eq!(
			earcut_no_holes(&[1.0, 0.0, 0.5, 0.5, 0.0, 1.0, 0.0, 0.0,]),
			[1, 2, 3, 1, 3, 0]
		);
	}

	#[test]
	fn big_circle() {
		let nv = 5000;
		let mut vertices = (0..nv)
			.map(|i| (i as f64) / (nv as f64) * f64::PI() * 1.5)
			.flat_map(|i| once(i.cos()).chain(once(i.sin())))
			.collect::<Vec<_>>();
		vertices.push(0.0);
		vertices.push(0.0);
		let triangles;
		{
			time_test!();
			triangles = earcut_no_holes(&vertices);
		}
		write_svg(&vertices, &triangles, "-1.2 -1.2 2.4 2.4").unwrap();
		assert_eq!(triangles.len(), 3 * (vertices.len() / 2 - 2));
	}
	#[test]
	fn cycle_windows_test() {
		let repeating = itertools::unfold(-1, |x| {
			*x = (*x + 1) % 4;
			Some(*x)
		});
		let mut actual = Vec::new();
		cycle_windows(
			repeating,
			|x| x,
			|arr| {
				actual.push(arr);
				None::<()>
			},
		);
		assert_eq!(actual, vec![[0, 1, 2], [1, 2, 3], [2, 3, 0], [3, 0, 1]]);
	}

	#[test]
	fn tt() {
		let w = 120;
		let h = 80;
		for y in 0..h {
			for x in 0..w {
				print!(
					"{}",
					if foo(
						lerp(-8.0, 8.0, x as f64 / w as f64),
						-lerp(-8.0, 8.0, y as f64 / h as f64),
						3.0,
						0.5
					) < 0.0
					{
						"X"
					} else {
						" "
					}
				);
			}
			println!("");
		}
	}

	use super::curve_point;
	use super::followAlgorithm2d;
	use super::write_point_svg;
	use itertools::Itertools;

	#[test]
	fn tt3() {
		let f = |x: f64, y: f64| foo(x, -y, 3.0, 0.5);
		let start_point = curve_point(&f, (3.0, -1.0));
		let (points, _) = followAlgorithm2d(&f, start_point, 0.1);
		// let (hole1_points, _) = followAlgorithm2d(&f, curve_point(&f, (0.1, 1.0)), 0.01);
		println!(
			"start_point: {:?} number_of_points  {}",
			start_point,
			points.len(),
		);
		// println!("{:?}", points);
		let mut vertices: Vec<f64> = Vec::new();
		for (x, y) in points {
			vertices.push(x);
			vertices.push(y);
		}
		write_point_svg(vertices.iter().cloned().tuples(), "-8 -8 16 16").unwrap();
		let triangles;
		{
			time_test!();
			triangles = earcut_no_holes(&vertices);
		}
		write_svg(&vertices, &triangles, "-8 -8 16 16").unwrap();
	}

	use super::archimedean;
	use super::pair_iter;
	use super::*;

	#[test]
	fn spiral() {
		// let points = archimedean(640, -4.0, f64::PI() * 8.0)
		let points = multi_archi(15.0, 3)
			// .map(|(x, y)| (-x, -y))
			// .chain(archimedean(640, f64::PI() * 8.0))
			.collect::<Vec<(f64, f64)>>();
		let points2 = points
			.iter()
			.cloned()
			.flat_map(pair_iter)
			.collect::<Vec<_>>();
		let view_box = "-30 -30 60 60 ";
		write_point_svg(points, view_box).unwrap();
		let triangles;
		{
			time_test!();
			triangles = earcut_no_holes(&points2);
		}
		write_svg(&points2, &triangles, view_box).unwrap();
	}

	#[test]
	fn menger() {
		let (points, _) = menger_wipe(3);
		println!("{:?}", points);
		write_point_svg(points, "-0.1 -0.1 1.2 1.2");
	}

	fn find_hole_bridge_test() {
		// find_hole_bridge(, hole_left_most: usize, e: impl Fn(usize) -> (f64, f64))
	}
}
fn lerp(a: f64, b: f64, t: f64) -> f64 {
	(1.0 - t) * a + t * b
}
use palette::{Hsl, Srgb};
use rand::Rng;
fn write_svg(vertices: &[f64], triangles: &[usize], view_box: &str) -> std::io::Result<()> {
	// return Ok(());
	let mut file = File::create("foo.svg")?;
	writeln!(
		file,
		"<?xml version='1.0' encoding='UTF-8' standalone='no'?>"
	)?;
	writeln!(
		file,
		"<svg height='1000' width='1000'  viewBox='{}' xmlns='http://www.w3.org/2000/svg'>",
		view_box
	)?;
	for (v, (i, j, k)) in triangles.iter().tuples().enumerate() {
		let color = Hsl::new(v as f64 / (triangles.len() as f64 / 3.0) * 270.0, 1.0, 0.5);
		let c2 = Srgb::from(color);
		writeln!(
			file,
			"<polygon points='{},{} {},{} {},{}' style='fill:rgb({}%,{}%,{}%);stroke:black;stroke-width:0.01;' />",
			vertices[i * 2],
			vertices[i * 2 + 1],
			vertices[j * 2],
			vertices[j * 2 + 1],
			vertices[k * 2],
			vertices[k * 2 + 1],
			c2.red * 100.0,
			c2.green * 100.0,
			c2.blue * 100.0,
		)?;
	}
	for i in 0..vertices.len() / 2 {
		writeln!(
			file,
			"<text x='{}' y='{}' font-size='0.1'>{}</text>",
			vertices[i * 2],
			vertices[i * 2 + 1],
			i,
		)?;
	}
	writeln!(file, "</svg>")?;
	Ok(())
}
fn write_point_svg(
	vertices: impl IntoIterator<Item = (f64, f64)>,
	view_box: &str,
) -> std::io::Result<()> {
	// return Ok(());
	let mut file = File::create("points.svg")?;
	writeln!(
		file,
		"<?xml version='1.0' encoding='UTF-8' standalone='no'?>"
	)?;
	writeln!(
		file,
		"<svg height='1000' width='1000'  viewBox='{}' xmlns='http://www.w3.org/2000/svg'>",
		view_box
	)?;
	writeln!(
		file,
		"<polygon stroke-width='0.01' stroke='black' fill='none' points='"
	)?;
	for (a, b) in vertices {
		writeln!(file, "{},{}", a, b)?;
	}
	writeln!(file, "' /></svg>")?;
	Ok(())
}
