use bit_vec::BitVec;
use itertools::Itertools;
use num_traits::real::Real;
use ord_subset::OrdSubsetIterExt;
use ord_subset::OrdSubsetSliceExt;
use std::fs::File;
use std::io::Write;
use std::iter::FromIterator;

#[must_use]
pub fn earcut(data: &[P2], hole_indices: &[usize]) -> (Vec<[usize; 3]>, Vec<P2>, Vec<P2>) {
	let nv = data.len();
	println!("nv {}", nv);
	let mut e = Vec::with_capacity(nv + hole_indices.len() * 2);
	e.extend_from_slice(data);
	println!("data2.len() {}", e.len());

	let hole_ranges = || {
		hole_indices.iter().enumerate().map(|(i, &start)| {
			let end = hole_indices.get(i + 1).cloned().unwrap_or(data.len());
			(start, end, end - start)
		})
	};
	let outer_loop_end = hole_indices.get(0).cloned().unwrap_or(data.len());

	let mut next_unused_index = Vec::with_capacity(e.capacity());
	for i in 0..nv {
		next_unused_index.push(i + 1);
	}
	for (start, end, _) in once((0, outer_loop_end, outer_loop_end)).chain(hole_ranges()) {
		next_unused_index[end - 1] = start;
	}

	// struct TriangleCycleIterator<'a, T: Copy> {
	// 	w: (T, T, T),
	// 	i: std::slice::Iter<'a, T>,
	// }
	// impl<'a, T: Copy> TriangleCycleIterator<'a, T> {
	// 	fn new(data: &[T]) -> TriangleCycleIterator<T> {
	// 		TriangleCycleIterator {
	// 			w: (
	// 				data[data.len() - 3],
	// 				data[data.len() - 2],
	// 				data[data.len() - 1],
	// 			),
	// 			i: data.iter(),
	// 		}
	// 	}
	// }
	// impl<'a, T: Copy> Iterator for TriangleCycleIterator<'a, T> {
	// 	type Item = (T, T, T);
	// 	fn next(&mut self) -> Option<(T, T, T)> {
	// 		self.i.next().map(|next| {
	// 			self.w[0] = self.w[0];
	// 			self.w[0] = self.w.2;
	// 			self.w.2 = *next;

	// 			self.w
	// 		})
	// 	}
	// }

	let mut is_convex = BitVec::with_capacity(e.capacity());
	for (start, _, d) in once((0, outer_loop_end, outer_loop_end)).chain(hole_ranges()) {
		for i in 0..d {
			is_convex.push(triangle_is_convex(
				e[start + (i + d - 1) % d],
				e[start + i],
				e[start + (i + 1) % d],
			))
		}
	}

	// (hole_range, hole_left_most_i, hole_left_most_x)
	let mut hole_left_mosts = hole_ranges()
		.map(|range| {
			let (start, end, _) = range;
			(start..end)
				.map(|i| (range, i, data[i][0]))
				.ord_subset_min_by_key(|&(_, _, x)| x)
				.unwrap()
		})
		.collect::<Vec<_>>();
	hole_left_mosts.ord_subset_sort_by_key(|&(_, _, x)| x);

	for ((start, _, length), hole_i, _) in hole_left_mosts {
		let hole_i_prev = start + ((hole_i - start) + 1) % length;
		let (bridge_i_prev, bridge_i) =
			find_hole_bridge(&next_unused_index, &is_convex, hole_i, |i| e[i]).unwrap();

		// bridge_i_prev -> bridge_i -> bridge_i_next and
		// hole_i_prev -> hole_i -> hole_i_next
		// become
		// bridge_i_prev -> bridge_i -> hole_i2 -> hole_i_next
		// hole_i_prev -> hole_i -> bridge_i2 -> bridge_i_next
		let bridge_i_next = next_unused_index[bridge_i];
		let hole_i_next = next_unused_index[hole_i];

		let hole_i2 = e.len();
		next_unused_index[bridge_i] = hole_i2;
		e.push(e[hole_i]);
		next_unused_index.push(hole_i_next);

		let bridge_i2 = e.len();
		next_unused_index[hole_i] = bridge_i2;
		e.push(e[bridge_i]);
		next_unused_index.push(bridge_i_next);

		println!(
			"bridge_i_prev {} bi {} bi_next {}",
			bridge_i_prev, bridge_i, hole_i2
		);

		is_convex.set(
			bridge_i,
			triangle_is_convex(e[bridge_i_prev], e[bridge_i], e[hole_i2]),
		);
		is_convex.set(
			hole_i,
			triangle_is_convex(e[hole_i_prev], e[hole_i], e[bridge_i2]),
		);
		is_convex.push(triangle_is_convex(e[bridge_i], e[hole_i2], e[hole_i_next]));
		is_convex.push(triangle_is_convex(
			e[hole_i],
			e[bridge_i2],
			e[bridge_i_next],
		));
	}

	let mut data3 = Vec::new();
	let mut i = 0;
	loop {
		let p = e[i];
		data3.push(p);
		i = next_unused_index[i];
		if i == 0 {
			break;
		}
	}
	println!();
	println!();
	println!("data3 {:?}", data3);
	let (triangles, left_over) = earcut_no_holes(&data3);
	(triangles, data3, left_over)
}

fn archimedean(steps: usize, start: f64, turns: f64) -> impl Iterator<Item = P2> {
	assert_ne!(steps, 0, "Steps should be at least 1.");
	let turns_per_step = (turns - start) / steps as f64;
	(0..=steps).map(move |i| {
		let rad = start + (i as f64 * turns_per_step);
		[rad * rad.cos(), rad * rad.abs().sin()]
	})
}

static PI: f64 = std::f64::consts::PI;
static TAU: f64 = PI * 2.0;

type M2 = [f64; 4];
fn apply(p: P2, m: M2) -> P2 {
	let [x, y] = p;
	let [a, b, c, d] = m;
	[a * x + b * y, c * x + d * y]
}

fn multi_archi(center_distance: f64, num: usize) -> impl Iterator<Item = P2> {
	// const center_distance: f64 = 14.0;
	// const num: usize = 5;
	fn rot(x: f64) -> M2 {
		[x.cos(), -x.sin(), x.sin(), x.cos()]
	}

	let closures = (0..num).map(move |i| {
		move |[x, y]: P2| {
			let trans = [x + center_distance, y];
			let mat = rot(TAU * i as f64 / num as f64);
			apply(trans, mat)
		}
	});
	closures.flat_map(|c| archimedean(128, -7.0, 9.0).map(c))
}

fn earcut_no_holes(data: &[P2]) -> (Vec<[usize; 3]>, Vec<P2>) {
	let modulo = |a: usize, b: usize| -> usize { ((a % b) + b) % b };
	let nv = data.len();
	let e = |i: usize| data[i];
	let mut triangles = Vec::with_capacity(3 * (nv - 2));
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
	for &[x, y] in data.iter() {
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

	let z_order_val = |[x, y]: P2| {
		interleave(
			(u16::max_value() as f64 * (x - min_x) * width_inv) as u16,
			(u16::max_value() as f64 * (y - min_y) * height_inv) as u16,
		)
	};
	let mut z_order_values = Vec::from_iter(data.iter().cloned().map(z_order_val).enumerate());
	z_order_values.sort_by_key(|(_, z): &(_, u32)| *z);
	let z_index_arr = {
		let mut x = vec![0; nv];
		for (i, j) in z_order_values.iter().map(|(j, _)| *j).enumerate() {
			x[j] = i;
		}
		x
	};

	fn is_ear(
		is_convex: &BitVec,
		is_used: &BitVec,
		e: impl Fn(usize) -> P2,
		i: usize,
		j: usize,
		k: usize,
	) -> bool {
		println!("is_ear ijk {} {} {} is_convex {:?}", i, j, k, is_convex[j]);
		let a = e(i);
		let b = e(j);
		let c = e(k);
		println!("a={:?} b={:?} c={:?}", a, b, c);
		is_convex[j]
			&& is_used.iter().enumerate().all(|(vi, used)| {
				let val = used
					|| vi == i || vi == j
					|| vi == k || !triangle_contains_point(a, b, c, e(vi))
					|| is_convex[vi];
				if !val {
					println!(
						"vi = {} e(vi) = {:?} is_convex[vi] = {}",
						vi,
						e(vi),
						is_convex[vi],
					);
				};
				val
			})
	};

	fn is_ear_zorder(
		is_convex: &BitVec,
		is_used: &BitVec,
		z_order_values: &[(usize, u32)],
		z_order_val: impl Fn(P2) -> u32,
		z_index_arr: &[usize],
		e: impl Fn(usize) -> P2,
		i: usize,
		j: usize,
		k: usize,
	) -> bool {
		if !is_convex[j] {
			return false;
		}

		// println!("is_ear_zorder {} {} {}", i, j, k);

		let nv = z_order_values.len();

		let [ax, ay] = e(i);
		let [bx, by] = e(j);
		let [cx, cy] = e(k);

		let min_triang_x = ax.min(bx).min(cx);
		let max_triang_x = ax.max(bx).max(cx);
		let min_triang_y = ay.min(by).min(cy);
		let max_triang_y = ay.max(by).max(cy);

		let min_triang_z = z_order_val([min_triang_x, min_triang_y]);
		let max_triang_z = z_order_val([max_triang_x, max_triang_y]);
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
				|| !triangle_contains_point([ax, ay], [bx, by], [cx, cy], e(vi))
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

	let mut left_over = Vec::new();

	while vertices_left >= 3 {
		{
			let [_, i, j, k] = current_indexes;
			// println!("{:?}", current_indexes);
			let start_end_index = k;
			// find next ear
			// while !is_ear_zorder(
			while !is_ear(
				&is_convex,
				&is_used,
				// &z_order_values,
				// z_order_val,
				// &z_index_arr,
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
					let start = is_used
						.iter()
						.enumerate()
						.filter(|(i, used)| !used)
						.next()
						.unwrap()
						.0;
					let mut i = start;
					loop {
						left_over.push(e(i));
						println!("{:?}", (i, e(i)));
						i = next_unused_index[i];
						if i == start {
							break;
						}
					}
					return (triangles, left_over);
				}
			}
		}
		{
			let [_, i, j, k] = current_indexes;
			triangles.push([i, j, k]);
			is_used.set(j, true);
			vertices_left = vertices_left - 1;
			println!("added, left: {}", vertices_left);
			println!();
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
		current_indexes[0] = current_indexes[1];
		current_indexes[1] = current_indexes[2];
		current_indexes[2] = current_indexes[3];
		current_indexes[3] = unused_index_after(&next_unused_index, current_indexes[2]);
	}
	(triangles, left_over)
}

fn double_area([ax, ay]: P2, [bx, by]: P2, [cx, cy]: P2) -> f64 {
	(bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
}

fn triangle_is_convex(a: P2, b: P2, c: P2) -> bool {
	// println!("{} {} {} {} {} {} ", ax, ay, bx, by, cx, cy);
	// println!("{}", (bx - ax) * (cy - ay) - (by - ay) * (cx - ax));
	double_area(a, b, c) > 0.0
}

fn triangle_contains_point([ax, ay]: P2, [bx, by]: P2, [cx, cy]: P2, [px, py]: P2) -> bool {
	let result = (cx - px) * (ay - py) - (ax - px) * (cy - py) >= 0.0
		&& (ax - px) * (by - py) - (bx - px) * (ay - py) >= 0.0
		&& (bx - px) * (cy - py) - (cx - px) * (by - py) >= 0.0;
	// println!(
	// 	"triangle_contains_point {},{} {},{} {},{}  {} {} => {}",
	// 	ax, ay, bx, by, cx, cy, px, py, result
	// );
	result
}

fn foo(x: f64, y: f64, a: f64, c: f64) -> f64 {
	x.powi(2) + y.powi(2)
		- a * ((x - y).sin().powi(3).exp() + (-x - y).sin().powi(3).exp()).powi(2)
		- c
}

fn ddd(f: impl Fn(f64, f64) -> f64, v: P2) -> P2 {
	let eps = 1e-6;
	let [x, y] = v;
	let fxy = f(x, y);
	[(f(x + eps, y) - fxy) / eps, (f(x, y + eps) - fxy) / eps]
}

fn curve_point(f: &impl Fn(f64, f64) -> f64, start: P2) -> P2 {
	let mut p = start;
	for _ in 0..8 {
		let [px, py] = p;
		let fp = f(px, py);
		let [dfpdx, dfpdy] = ddd(f, p);
		let scale = fp / (dfpdx * dfpdx + dfpdy * dfpdy);
		p = [px - scale * dfpdx, py - scale * dfpdy];
		// println!("p {:?}       f {:?}", p, f[px, py]);
	}
	p
}

fn distance(a: P2, b: P2) -> f64 {
	let [ax, ay] = a;
	let [bx, by] = b;
	length([ax - bx, ay - by])
}

type P2 = [f64; 2];

fn to_length(v: P2, nl: f64) -> P2 {
	let l = length(v);
	let [x, y] = v;
	[x * nl / l, y * nl / l]
}

fn length(v: P2) -> f64 {
	let [x, y] = v;
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
			min: [INF, INF],
			max: [-INF, -INF],
		}
	}
	fn add_iter(self: &mut Self, iter: impl IntoIterator<Item = P2>) {
		for [x, y] in iter {
			self.min[0] = self.min[0].min(x);
			self.max[0] = self.max[0].max(x);

			self.min[0] = self.min[0].min(y);
			self.max[0] = self.max[0].max(y);
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
	e: impl Fn(usize) -> P2,
) -> Option<(usize, usize)> {
	let [hole_point_x, hole_point_y] = e(hole_left_most); // hole x/y
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
		|[(i_prev, _), (i, [px, py]), (i_next, [p_next_x, p_next_y])]| {
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

		let (stop, _) = m;
		let mut tan_min = INF;

		// p = m.next;
		let (mut mi_prev, mut mi) = m;
		let (mut pi_prev, mut pi) = m;
		let [mx, my] = e(mi);

		loop {
			pi_prev = pi;
			pi = next_unused_index[pi];

			let [px, py] = e(pi);
			if hole_point_x >= px
				&& px >= mx && hole_point_x != px
				&& triangle_contains_point(
					[
						if hole_point_y < my { hole_point_x } else { qx },
						hole_point_y,
					],
					[mx, my],
					[
						if hole_point_y < my { qx } else { hole_point_x },
						hole_point_y,
					],
					[px, py],
				) {
				let tan = (hole_point_y - py).abs() / (hole_point_x - px); // tangential

				if (tan < tan_min || (tan == tan_min && px > mx))
					&& locally_inside(
						e(pi_prev),
						e(pi),
						e(next_unused_index[pi]),
						[hole_point_x, hole_point_y],
					) {
					mi = pi;
					mi_prev = pi_prev;
					tan_min = tan;
				}
			}

			if pi == stop {
				break;
			}
		}

		(mi_prev, mi)
	})
}

fn mm([a, b, c, d]: M2, [e, f, g, h]: M2) -> M2 {
	[a * e + b * g, a * f + b * h, c * e + d * g, c * f + d * h]
}

fn menger_wipe(iterations: u32) -> (Vec<P2>, Vec<usize>) {
	fn fn_menger_recursive(
		points: &mut Vec<P2>,
		hole_starts: &mut Vec<usize>,
		lvl: u32,
		[ox, oy]: P2,
		s: f64,
	) {
		hole_starts.push(points.len());
		for [x, y] in &[
			[1.0 / 3.0, 1.0 / 3.0],
			[1.0 / 3.0, 2.0 / 3.0],
			[2.0 / 3.0, 2.0 / 3.0],
			[2.0 / 3.0, 1.0 / 3.0],
		] {
			points.push([ox + s * x, oy + s * y]);
		}

		if 0 != lvl {
			for i in 0..3 {
				for j in 0..3 {
					if j != 1 || i != 1 {
						if lvl == 3 && i == 0 && j == 0
							|| lvl == 2 && j <= 1 && !(i == 0 && j == 1)
							|| lvl != 3 && lvl != 2
						{
							fn_menger_recursive(
								points,
								hole_starts,
								lvl - 1,
								[ox + (i as f64 / 3.0) * s, oy + (j as f64 / 3.0) * s],
								s / 3.0,
							);
						}
					}
				}
			}
		}
	};

	let hole_count = (8usize.pow(1 + iterations) - 1) / 7;
	let mut points = Vec::with_capacity((1 + hole_count) * 4);
	let mut hole_starts = Vec::with_capacity(hole_count);
	points.push([0.0, 0.0]);
	points.push([1.0, 0.0]);
	points.push([1.0, 1.0]);
	points.push([0.0, 1.0]);
	fn_menger_recursive(&mut points, &mut hole_starts, iterations, [0.0, 0.0], 1.0);
	(points, hole_starts)
}

// check if k is inside the angle (i.e. to the left of the angle spanned by a -> b -> c)
fn locally_inside(a: P2, b: P2, c: P2, k: P2) -> bool {
	// if a b c form a concave angle (on the left), then k must be left of both the lines spanned
	// by a b and b c for it to be inside the angle.
	// otherwise, is suffices for it to be left of either.
	let mine = if double_area(a, b, c) > 0.0 {
		double_area(a, b, k) >= 0.0 && double_area(b, c, k) >= 0.0
	} else {
		double_area(a, b, k) >= 0.0 || double_area(b, c, k) >= 0.0
	};
	mine
}

fn followAlgorithm2d(
	ic: &impl Fn(f64, f64) -> f64,
	start_p: P2,
	step_length: f64,
) -> (Vec<P2>, Vec<P2>) {
	let [ddx, ddy] = ddd(ic, start_p);
	let start_tangent = to_length([-ddy, ddx], step_length);
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
		let [px, py] = p;
		let [tx, ty] = tangent;
		let search_start = [px + tx, py + ty];
		let new_p = curve_point(&ic, search_start);
		let [ddx, ddy] = ddd(ic, new_p);
		let new_tangent = to_length([-ddy, ddx], step_length);
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

	use super::*;

	use num_traits::float::FloatConst;
	use std::iter::once;

	#[test]
	fn _3vertices() {
		// assert_eq!(earcut_no_holes(&[0.0, 0.0, 1.0, 0.0, 0.0, 1.0,]), [1, 2, 0]);
	}

	#[test]
	fn _4vertices() {
		// assert_eq!(
		// 	earcut_no_holes(&[0.0, 0.0, 1.0, 0.0, 0.5, 0.5, 0.0, 1.0,]),
		// 	[2, 3, 0, 2, 0, 1]
		// );
		// assert_eq!(
		// 	earcut_no_holes(&[1.0, 0.0, 0.5, 0.5, 0.0, 1.0, 0.0, 0.0,]),
		// 	[1, 2, 3, 1, 3, 0]
		// );
	}

	#[test]
	fn big_circle() {
		let nv = 5000;
		let vertices = (0..nv)
			.map(|i| (i as f64) / (nv as f64) * f64::PI() * 1.5)
			.map(|i| [i.cos(), i.sin()])
			.chain(once([0.0, 0.0]))
			.collect::<Vec<_>>();
		let triangles;
		{
			time_test!();
			triangles = earcut_no_holes(&vertices).0;
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
						0.5,
					) < 0.0
					{
						"X"
					} else {
						" "
					}
				);
			}
			println!();
		}
	}

	#[test]
	fn tt3() {
		let f = |x: f64, y: f64| foo(x, -y, 3.0, 0.5);
		let start_point = curve_point(&f, [3.0, -1.0]);
		let (points, _) = followAlgorithm2d(&f, start_point, 0.1);
		// let (hole1_points, _) = followAlgorithm2d(&f, curve_point(&f, [0.1, 1.0]), 0.01);
		println!(
			"start_point: {:?} number_of_points  {}",
			start_point,
			points.len(),
		);
		// println!("{:?}", points);
		// write_point_svg(vertices.iter().cloned(), "-8 -8 16 16").unwrap();
		let triangles;
		{
			time_test!();
			triangles = earcut_no_holes(&points).0;
		}
		write_svg(&points, &triangles, "-8 -8 16 16").unwrap();
	}

	#[test]
	fn spiral() {
		// let points = archimedean(640, -4.0, f64::PI() * 8.0)
		let points = multi_archi(15.0, 3)
			// .map(|[x, y]| (-x, -y))
			// .chain(archimedean(640, f64::PI() * 8.0))
			.collect::<Vec<P2>>();
		let view_box = "-30 -30 60 60 ";
		// write_point_svg(points, view_box).unwrap();
		let triangles;
		{
			time_test!();
			triangles = earcut_no_holes(&points).0;
		}
		write_svg(&points, &triangles, view_box).unwrap();
	}

	#[test]
	fn triangle_contains_point_test() {
		assert_eq!(
			true,
			triangle_contains_point([0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 0.0])
		);
	}

	#[test]
	fn menger1() {
		let (points, hole_offsets) = menger_wipe(3);
		for p in points {
			let [x, y] = p;
			if x < 1.0 / 9.0 + 1e-6 && y > 1.0 / 9.0 - 1e-6 && y < 2.0 / 9.0 + 1e-6 {
				println!("{:?}", p);
			}
		}
		return;
		// println!("{:?}", points);
		let (triangles, data3, left_over) = earcut(
			&points
				.iter()
				.map(|[x, y]| [(x * 720.0), (y * 720.0)])
				.collect::<Vec<_>>(),
			&hole_offsets,
		);
		let view_box = "-1 -1 730 730";
		write_point_svg(&data3, view_box).unwrap();
		write_point_svg(&left_over, view_box).unwrap();
		// let (t2, hlo2) = earcut_no_holes(&data5);
		write_svg(&data3, &triangles, view_box).unwrap();
		// write_svg(&hlo2, &t2, view_box);
		// assert_eq!(24, triangles.len());
	}

	#[test]
	fn menger2() {
		let (points, hole_offsets) = menger_wipe(3);
		// println!("{:?}", points);
		let points_moved = points
			.iter()
			.map(|[x, y]| [x * 1.0 / 3.0 / 3.0, 1.0 / 9.0 + y * 1.0 / 3.0 / 3.0])
			.collect::<Vec<_>>();
		let points_scaled = points
			.iter()
			.map(|[x, y]| [(x * 720.0), (y * 720.0)])
			.collect::<Vec<_>>();
		// return;
		let (triangles, data3, left_over) = earcut(&points_scaled, &hole_offsets);
		let x: [usize; 6] = [74, 70, 122, 123, 130, 131];
		for &i in x.iter() {
			println!("{:?},", data3[i]);
		}
		let view_box = "-1 -1 730 730";
		write_point_svg(&data3, view_box).unwrap();
		write_point_svg(&left_over, view_box).unwrap();
		// let (t2, hlo2) = earcut_no_holes(&data5);
		write_svg(&data3, &triangles, view_box).unwrap();
		// write_svg(&hlo2, &t2, view_box);
		// assert_eq!(24, triangles.len());
	}

	#[test]
	fn menger3() {
		let p = vec![[204.44444444444443, 71.11111111111111]];
		let (points, hole_offsets) = menger_wipe(3);
		// println!("{:?}", points);
		let points_moved = points
			.iter()
			.map(|[x, y]| [x * 1.0 / 3.0 / 3.0, 1.0 / 9.0 + y * 1.0 / 3.0 / 3.0])
			.collect::<Vec<_>>();
		let points_scaled = points
			.iter()
			.map(|[x, y]| [(x * 720.0), (y * 720.0)])
			.collect::<Vec<_>>();
		for p in points_moved {
			println!("{:?}", p);
		}
		// return;
		let (triangles, data3, left_over) = earcut(&points_scaled, &hole_offsets);
		let view_box = "-1 -1 730 730";
		write_point_svg(&data3, view_box).unwrap();
		// write_point_svg(&left_over, view_box).unwrap();
		// let (t2, hlo2) = earcut_no_holes(&data5);
		write_svg(&data3, &triangles, view_box).unwrap();
		// write_svg(&hlo2, &t2, view_box);
		// assert_eq!(24, triangles.len());
	}

	fn find_hole_bridge_test() {
		// find_hole_bridge(, hole_left_most: usize, e: impl Fn(usize) -> P2)
	}

	#[test]
	fn tt4() {
		let p = vec![
			[123.25925925925925, 90.07407407407408],
			[85.33333333333333, 85.33333333333333],
			[85.33333333333333, 42.666666666666664],
			[90.07407407407408, 47.407407407407405],
			[90.07407407407408, 52.148148148148145],
			[99.55555555555554, 56.888888888888886],
			[99.55555555555554, 71.11111111111111],
			[113.77777777777777, 71.11111111111111],
			[113.77777777777777, 56.888888888888886],
			[118.5185185185185, 61.629629629629626],
			[118.5185185185185, 66.37037037037037],
			[123.25925925925925, 66.37037037037037],
			[123.25925925925925, 61.629629629629626],
			[113.77777777777777, 56.888888888888886],
			[99.55555555555554, 56.888888888888886],
			[90.07407407407408, 52.148148148148145],
			[94.81481481481481, 52.148148148148145],
			[94.81481481481481, 47.407407407407405],
			[104.29629629629628, 47.407407407407405],
			[104.29629629629628, 52.148148148148145],
			[109.03703703703702, 52.148148148148145],
			[109.03703703703702, 47.407407407407405],
			[118.5185185185185, 47.407407407407405],
			[118.5185185185185, 52.148148148148145],
			[123.25925925925925, 52.148148148148145],
			[123.25925925925925, 47.407407407407405],
			[118.5185185185185, 47.407407407407405],
			[109.03703703703702, 47.407407407407405],
			[104.29629629629628, 47.407407407407405],
			[94.81481481481481, 47.407407407407405],
			[90.07407407407408, 47.407407407407405],
			[85.33333333333333, 42.666666666666664],
			[90.07407407407408, 75.85185185185185],
			[90.07407407407408, 80.5925925925926],
			[94.81481481481481, 80.5925925925926],
			[94.81481481481481, 75.85185185185185],
			[104.29629629629628, 75.85185185185185],
			[104.29629629629628, 80.5925925925926],
			[109.03703703703702, 80.5925925925926],
			[109.03703703703702, 75.85185185185185],
			[118.5185185185185, 75.85185185185185],
			[118.5185185185185, 80.5925925925926],
			[123.25925925925925, 80.5925925925926],
			[123.25925925925925, 75.85185185185185],
			[118.5185185185185, 75.85185185185185],
			[109.03703703703702, 75.85185185185185],
			[104.29629629629628, 75.85185185185185],
			[94.81481481481481, 75.85185185185185],
			[90.07407407407408, 75.85185185185185],
			[85.33333333333333, 42.666666666666664],
			[90.07407407407408, 61.629629629629626],
			[90.07407407407408, 66.37037037037037],
			[94.81481481481481, 66.37037037037037],
			[94.81481481481481, 61.629629629629626],
			[90.07407407407408, 61.629629629629626],
			[85.33333333333333, 42.666666666666664],
			[42.666666666666664, 42.666666666666664],
			[123.25925925925925, 37.925925925925924],
			[128.0, 128.0],
		];
		let view_box = "-1 -1 130 130";
		write_point_svg(&p, view_box).unwrap();
		let triangles;
		{
			time_test!();
			triangles = earcut_no_holes(&p).0;
		}
		write_svg(&p, &triangles, view_box).unwrap();
	}
	#[test]
	fn tt6() {
		let p = vec![
			[53.33333333333333, 53.33333333333333],
			[44.44444444444444, 62.22222222222222],
			[17.77777777777778, 62.22222222222222],
			[8.88888888888889, 44.44444444444444],
			[26.666666666666664, 53.33333333333333],
		];
		let view_box = "-1 -1 730 730";
		write_point_svg(&p, view_box).unwrap();
		let triangles;
		{
			time_test!();
			triangles = earcut_no_holes(&p).0;
		}
		write_svg(&p, &triangles, view_box).unwrap();
	}

	#[test]
	fn tt5() {
		let p = vec![
			[711.1111111111111, 222.2222222222222],
			[693.3333333333333, 213.33333333333331],
			[693.3333333333333, 186.66666666666666],
			[702.2222222222222, 195.55555555555554],
			[702.2222222222222, 204.44444444444443],
			[711.1111111111111, 204.44444444444443],
			[711.1111111111111, 195.55555555555554],
			[693.3333333333333, 186.66666666666666],
			[151.11111111111111, 177.77777777777777],
			[151.11111111111111, 168.88888888888889],
			[168.88888888888889, 168.88888888888889],
			[168.88888888888889, 177.77777777777777],
			[177.77777777777777, 177.77777777777777],
			[177.77777777777777, 168.88888888888889],
			[195.55555555555554, 168.88888888888889],
			[195.55555555555554, 177.77777777777777],
			[204.44444444444443, 177.77777777777777],
			[204.44444444444443, 168.88888888888889],
			[222.2222222222222, 168.88888888888889],
			[222.2222222222222, 177.77777777777777],
			[231.1111111111111, 177.77777777777777],
			[231.1111111111111, 168.88888888888889],
			[248.88888888888889, 168.88888888888889],
			[248.88888888888889, 177.77777777777777],
			[257.77777777777777, 177.77777777777777],
			[257.77777777777777, 168.88888888888889],
			[275.55555555555554, 168.88888888888889],
			[275.55555555555554, 177.77777777777777],
			[284.44444444444446, 177.77777777777777],
			[284.44444444444446, 168.88888888888889],
			[302.22222222222223, 168.88888888888889],
			[302.22222222222223, 177.77777777777777],
			[311.1111111111111, 177.77777777777777],
			[311.1111111111111, 168.88888888888889],
			[328.88888888888886, 168.88888888888889],
			[328.88888888888886, 177.77777777777777],
			[337.77777777777777, 177.77777777777777],
			[337.77777777777777, 168.88888888888889],
			[355.55555555555554, 168.88888888888889],
			[355.55555555555554, 177.77777777777777],
			[364.44444444444446, 177.77777777777777],
			[364.44444444444446, 168.88888888888889],
			[382.2222222222222, 168.88888888888889],
			[382.2222222222222, 177.77777777777777],
			[391.1111111111111, 177.77777777777777],
			[391.1111111111111, 168.88888888888889],
			[408.88888888888897, 168.88888888888889],
			[408.88888888888897, 177.77777777777777],
			[417.77777777777777, 177.77777777777777],
			[417.77777777777777, 168.88888888888889],
			[435.55555555555554, 168.88888888888889],
			[435.55555555555554, 177.77777777777777],
			[444.4444444444444, 177.77777777777777],
			[444.4444444444444, 168.88888888888889],
			[462.2222222222222, 168.88888888888889],
			[462.2222222222222, 177.77777777777777],
			[471.11111111111114, 177.77777777777777],
			[471.11111111111114, 168.88888888888889],
			[488.8888888888888, 168.88888888888889],
			[488.8888888888888, 177.77777777777777],
			[497.77777777777777, 177.77777777777777],
			[497.77777777777777, 168.88888888888889],
			[515.5555555555557, 168.88888888888889],
			[515.5555555555557, 177.77777777777777],
			[524.4444444444445, 177.77777777777777],
			[524.4444444444445, 168.88888888888889],
			[542.2222222222223, 168.88888888888889],
			[542.2222222222223, 177.77777777777777],
			[551.1111111111111, 177.77777777777777],
			[551.1111111111111, 168.88888888888889],
			[568.8888888888889, 168.88888888888889],
			[568.8888888888889, 177.77777777777777],
			[577.7777777777777, 177.77777777777777],
			[577.7777777777777, 168.88888888888889],
			[595.5555555555554, 168.88888888888889],
			[595.5555555555554, 177.77777777777777],
			[604.4444444444443, 177.77777777777777],
			[604.4444444444443, 168.88888888888889],
			[622.2222222222221, 168.88888888888889],
			[622.2222222222221, 177.77777777777777],
			[631.1111111111111, 177.77777777777777],
			[631.1111111111111, 168.88888888888889],
			[648.8888888888889, 168.88888888888889],
			[648.8888888888889, 177.77777777777777],
			[657.7777777777777, 177.77777777777777],
			[657.7777777777777, 168.88888888888889],
			[675.5555555555555, 168.88888888888889],
			[675.5555555555555, 177.77777777777777],
			[684.4444444444443, 177.77777777777777],
			[684.4444444444443, 168.88888888888889],
			[702.2222222222222, 168.88888888888889],
			[702.2222222222222, 177.77777777777777],
			[711.1111111111111, 177.77777777777777],
			[711.1111111111111, 168.88888888888889],
			[702.2222222222222, 168.88888888888889],
			[684.4444444444443, 168.88888888888889],
			[675.5555555555555, 168.88888888888889],
			[657.7777777777777, 168.88888888888889],
			[648.8888888888889, 168.88888888888889],
			[631.1111111111111, 168.88888888888889],
			[622.2222222222221, 168.88888888888889],
			[604.4444444444443, 168.88888888888889],
			[595.5555555555554, 168.88888888888889],
			[577.7777777777777, 168.88888888888889],
			[568.8888888888889, 168.88888888888889],
			[551.1111111111111, 168.88888888888889],
			[542.2222222222223, 168.88888888888889],
			[524.4444444444445, 168.88888888888889],
			[515.5555555555557, 168.88888888888889],
			[497.77777777777777, 168.88888888888889],
			[488.8888888888888, 168.88888888888889],
			[471.11111111111114, 168.88888888888889],
			[462.2222222222222, 168.88888888888889],
			[444.4444444444444, 168.88888888888889],
			[435.55555555555554, 168.88888888888889],
			[417.77777777777777, 168.88888888888889],
			[408.88888888888897, 168.88888888888889],
			[391.1111111111111, 168.88888888888889],
			[382.2222222222222, 168.88888888888889],
			[364.44444444444446, 168.88888888888889],
			[355.55555555555554, 168.88888888888889],
			[337.77777777777777, 168.88888888888889],
			[328.88888888888886, 168.88888888888889],
			[311.1111111111111, 168.88888888888889],
			[302.22222222222223, 168.88888888888889],
			[284.44444444444446, 168.88888888888889],
			[275.55555555555554, 168.88888888888889],
			[257.77777777777777, 168.88888888888889],
			[248.88888888888889, 168.88888888888889],
			[231.1111111111111, 168.88888888888889],
			[222.2222222222222, 168.88888888888889],
			[204.44444444444443, 168.88888888888889],
			[195.55555555555554, 168.88888888888889],
			[177.77777777777777, 168.88888888888889],
			[168.88888888888889, 168.88888888888889],
			[151.11111111111111, 168.88888888888889],
			[142.22222222222223, 168.88888888888889],
			[8.88888888888889, 151.11111111111111],
			[17.77777777777778, 151.11111111111111],
			[35.55555555555556, 151.11111111111111],
			[44.44444444444444, 151.11111111111111],
			[62.22222222222222, 151.11111111111111],
			[71.11111111111111, 151.11111111111111],
			[71.11111111111111, 142.22222222222223],
			[71.11111111111111, 124.44444444444444],
			[71.11111111111111, 115.55555555555554],
			[53.33333333333333, 106.66666666666666],
			[17.77777777777778, 88.88888888888889],
			[35.55555555555556, 97.77777777777777],
			[44.44444444444444, 97.77777777777777],
			[62.22222222222222, 97.77777777777777],
			[71.11111111111111, 97.77777777777777],
			[80.0, 160.0],
			[160.0, 160.0],
			[320.0, 160.0],
			[400.0, 160.0],
			[559.9999999999999, 160.0],
			[640.0, 160.0],
			[640.0, 80.0],
			[648.8888888888889, 142.22222222222223],
			[648.8888888888889, 151.11111111111111],
			[657.7777777777777, 151.11111111111111],
			[657.7777777777777, 142.22222222222223],
			[675.5555555555555, 142.22222222222223],
			[675.5555555555555, 151.11111111111111],
			[684.4444444444443, 151.11111111111111],
			[684.4444444444443, 142.22222222222223],
			[702.2222222222222, 142.22222222222223],
			[702.2222222222222, 151.11111111111111],
			[711.1111111111111, 151.11111111111111],
			[711.1111111111111, 142.22222222222223],
			[702.2222222222222, 142.22222222222223],
			[684.4444444444443, 142.22222222222223],
			[675.5555555555555, 142.22222222222223],
			[657.7777777777777, 142.22222222222223],
			[648.8888888888889, 142.22222222222223],
			[640.0, 80.0],
			[648.8888888888889, 115.55555555555554],
			[648.8888888888889, 124.44444444444444],
			[657.7777777777777, 124.44444444444444],
			[657.7777777777777, 115.55555555555554],
			[648.8888888888889, 115.55555555555554],
			[640.0, 80.0],
			[648.8888888888889, 88.88888888888889],
			[648.8888888888889, 97.77777777777777],
			[666.6666666666666, 133.33333333333331],
			[693.3333333333333, 133.33333333333331],
			[693.3333333333333, 106.66666666666666],
			[702.2222222222222, 115.55555555555554],
			[702.2222222222222, 124.44444444444444],
			[711.1111111111111, 124.44444444444444],
			[711.1111111111111, 115.55555555555554],
			[693.3333333333333, 106.66666666666666],
			[666.6666666666666, 106.66666666666666],
			[648.8888888888889, 97.77777777777777],
			[657.7777777777777, 97.77777777777777],
			[675.5555555555555, 97.77777777777777],
			[684.4444444444443, 97.77777777777777],
			[702.2222222222222, 97.77777777777777],
			[711.1111111111111, 97.77777777777777],
			[711.1111111111111, 88.88888888888889],
			[640.0, 80.0],
			[44.44444444444444, 71.11111111111111],
			[44.44444444444444, 62.22222222222222],
			[62.22222222222222, 62.22222222222222],
			[62.22222222222222, 71.11111111111111],
			[71.11111111111111, 71.11111111111111],
			[71.11111111111111, 62.22222222222222],
			[88.88888888888889, 62.22222222222222],
			[88.88888888888889, 71.11111111111111],
			[97.77777777777777, 71.11111111111111],
			[97.77777777777777, 62.22222222222222],
			[115.55555555555554, 62.22222222222222],
			[115.55555555555554, 71.11111111111111],
			[124.44444444444444, 71.11111111111111],
			[124.44444444444444, 62.22222222222222],
			[142.22222222222223, 62.22222222222222],
			[142.22222222222223, 71.11111111111111],
			[151.11111111111111, 71.11111111111111],
			[151.11111111111111, 62.22222222222222],
			[168.88888888888889, 62.22222222222222],
			[168.88888888888889, 71.11111111111111],
			[177.77777777777777, 71.11111111111111],
			[177.77777777777777, 62.22222222222222],
			[195.55555555555554, 62.22222222222222],
			[195.55555555555554, 71.11111111111111],
			[204.44444444444443, 71.11111111111111],
			[204.44444444444443, 62.22222222222222],
			[222.2222222222222, 62.22222222222222],
			[222.2222222222222, 71.11111111111111],
			[231.1111111111111, 71.11111111111111],
			[231.1111111111111, 62.22222222222222],
			[248.88888888888889, 62.22222222222222],
			[248.88888888888889, 71.11111111111111],
			[257.77777777777777, 71.11111111111111],
			[257.77777777777777, 62.22222222222222],
			[275.55555555555554, 62.22222222222222],
			[275.55555555555554, 71.11111111111111],
			[284.44444444444446, 71.11111111111111],
			[284.44444444444446, 62.22222222222222],
			[302.22222222222223, 62.22222222222222],
			[302.22222222222223, 71.11111111111111],
			[311.1111111111111, 71.11111111111111],
			[311.1111111111111, 62.22222222222222],
			[328.88888888888886, 62.22222222222222],
			[328.88888888888886, 71.11111111111111],
			[337.77777777777777, 71.11111111111111],
			[337.77777777777777, 62.22222222222222],
			[355.55555555555554, 62.22222222222222],
			[355.55555555555554, 71.11111111111111],
			[364.44444444444446, 71.11111111111111],
			[364.44444444444446, 62.22222222222222],
			[382.2222222222222, 62.22222222222222],
			[382.2222222222222, 71.11111111111111],
			[391.1111111111111, 71.11111111111111],
			[391.1111111111111, 62.22222222222222],
			[408.88888888888897, 62.22222222222222],
			[408.88888888888897, 71.11111111111111],
			[417.77777777777777, 71.11111111111111],
			[417.77777777777777, 62.22222222222222],
			[435.55555555555554, 62.22222222222222],
			[435.55555555555554, 71.11111111111111],
			[444.4444444444444, 71.11111111111111],
			[444.4444444444444, 62.22222222222222],
			[462.2222222222222, 62.22222222222222],
			[462.2222222222222, 71.11111111111111],
			[471.11111111111114, 71.11111111111111],
			[471.11111111111114, 62.22222222222222],
			[488.8888888888888, 62.22222222222222],
			[488.8888888888888, 71.11111111111111],
			[497.77777777777777, 71.11111111111111],
			[497.77777777777777, 62.22222222222222],
			[515.5555555555557, 62.22222222222222],
			[515.5555555555557, 71.11111111111111],
			[524.4444444444445, 71.11111111111111],
			[524.4444444444445, 62.22222222222222],
			[542.2222222222223, 62.22222222222222],
			[542.2222222222223, 71.11111111111111],
			[551.1111111111111, 71.11111111111111],
			[551.1111111111111, 62.22222222222222],
			[568.8888888888889, 62.22222222222222],
			[568.8888888888889, 71.11111111111111],
			[577.7777777777777, 71.11111111111111],
			[577.7777777777777, 62.22222222222222],
			[595.5555555555554, 62.22222222222222],
			[595.5555555555554, 71.11111111111111],
			[604.4444444444443, 71.11111111111111],
			[604.4444444444443, 62.22222222222222],
			[622.2222222222221, 62.22222222222222],
			[622.2222222222221, 71.11111111111111],
			[631.1111111111111, 71.11111111111111],
			[631.1111111111111, 62.22222222222222],
			[648.8888888888889, 62.22222222222222],
			[648.8888888888889, 71.11111111111111],
			[657.7777777777777, 71.11111111111111],
			[657.7777777777777, 62.22222222222222],
			[675.5555555555555, 62.22222222222222],
			[675.5555555555555, 71.11111111111111],
			[684.4444444444443, 71.11111111111111],
			[684.4444444444443, 62.22222222222222],
			[702.2222222222222, 62.22222222222222],
			[702.2222222222222, 71.11111111111111],
			[711.1111111111111, 71.11111111111111],
			[711.1111111111111, 62.22222222222222],
			[702.2222222222222, 62.22222222222222],
			[684.4444444444443, 62.22222222222222],
			[675.5555555555555, 62.22222222222222],
			[657.7777777777777, 62.22222222222222],
			[648.8888888888889, 62.22222222222222],
			[631.1111111111111, 62.22222222222222],
			[622.2222222222221, 62.22222222222222],
			[604.4444444444443, 62.22222222222222],
			[595.5555555555554, 62.22222222222222],
			[577.7777777777777, 62.22222222222222],
			[568.8888888888889, 62.22222222222222],
			[551.1111111111111, 62.22222222222222],
			[542.2222222222223, 62.22222222222222],
			[524.4444444444445, 62.22222222222222],
			[515.5555555555557, 62.22222222222222],
			[497.77777777777777, 62.22222222222222],
			[488.8888888888888, 62.22222222222222],
			[471.11111111111114, 62.22222222222222],
			[462.2222222222222, 62.22222222222222],
			[444.4444444444444, 62.22222222222222],
			[435.55555555555554, 62.22222222222222],
			[417.77777777777777, 62.22222222222222],
			[408.88888888888897, 62.22222222222222],
			[391.1111111111111, 62.22222222222222],
			[382.2222222222222, 62.22222222222222],
			[364.44444444444446, 62.22222222222222],
			[355.55555555555554, 62.22222222222222],
			[337.77777777777777, 62.22222222222222],
			[328.88888888888886, 62.22222222222222],
			[311.1111111111111, 62.22222222222222],
			[302.22222222222223, 62.22222222222222],
			[284.44444444444446, 62.22222222222222],
			[275.55555555555554, 62.22222222222222],
			[257.77777777777777, 62.22222222222222],
			[248.88888888888889, 62.22222222222222],
			[231.1111111111111, 62.22222222222222],
			[222.2222222222222, 62.22222222222222],
			[204.44444444444443, 62.22222222222222],
			[195.55555555555554, 62.22222222222222],
			[177.77777777777777, 62.22222222222222],
			[168.88888888888889, 62.22222222222222],
			[151.11111111111111, 62.22222222222222],
			[142.22222222222223, 62.22222222222222],
			[124.44444444444444, 62.22222222222222],
			[115.55555555555554, 62.22222222222222],
			[97.77777777777777, 62.22222222222222],
			[88.88888888888889, 62.22222222222222],
			[71.11111111111111, 62.22222222222222],
			[62.22222222222222, 62.22222222222222],
			[44.44444444444444, 62.22222222222222],
			[8.88888888888889, 44.44444444444444],
			[26.666666666666664, 53.33333333333333],
			[53.33333333333333, 53.33333333333333],
			[106.66666666666666, 53.33333333333333],
			[133.33333333333331, 53.33333333333333],
			[186.66666666666666, 53.33333333333333],
			[213.33333333333331, 53.33333333333333],
			[266.66666666666663, 53.33333333333333],
			[293.3333333333333, 53.33333333333333],
			[346.66666666666663, 53.33333333333333],
			[373.3333333333333, 53.33333333333333],
			[426.66666666666663, 53.33333333333333],
			[453.33333333333337, 53.33333333333333],
			[506.6666666666667, 53.33333333333333],
			[533.3333333333333, 53.33333333333333],
			[586.6666666666665, 53.33333333333333],
			[613.3333333333333, 53.33333333333333],
			[666.6666666666666, 53.33333333333333],
			[693.3333333333333, 53.33333333333333],
			[711.1111111111111, 44.44444444444444],
			[720.0, 720.0],
		];
		let points2 = p
			.iter()
			.map(|&[x, y]| [x.round(), y.round()])
			.collect::<Vec<_>>();
		let view_box = "-1 -1 730 730";
		write_point_svg(&p, view_box).unwrap();
		let triangles;
		{
			time_test!();
			triangles = earcut_no_holes(&points2).0;
		}
		write_svg(&p, &triangles, view_box).unwrap();
	}
}

fn lerp(a: f64, b: f64, t: f64) -> f64 {
	(1.0 - t) * a + t * b
}

use palette::{Hsl, Srgb};
use rand::Rng;

fn write_svg(vertices: &[P2], triangles: &[[usize; 3]], view_box: &str) -> std::io::Result<()> {
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
	for (v, &[i, j, k]) in triangles.iter().enumerate() {
		let color = Hsl::new(v as f64 / (triangles.len() as f64 / 3.0) * 270.0, 1.0, 0.5);
		let c2 = Srgb::from(color);
		writeln!(
            file,
            "<polygon points='{},{} {},{} {},{}' style='fill:rgb({}%,{}%,{}%);stroke:black;stroke-width:0.01;' />",
            vertices[i][0],vertices[i][1],
            vertices[j][0],vertices[j][1],
            vertices[k][0],vertices[k][1],
            c2.red * 100.0,
            c2.green * 100.0,
            c2.blue * 100.0,
        )?;
	}
	for i in 0..vertices.len() {
		writeln!(
			file,
			"<text x='{}' y='{}' font-size='0.1'>{}</text>",
			vertices[i][0], vertices[i][1], i,
		)?;
	}
	writeln!(file, "</svg>")?;
	Ok(())
}

fn write_point_svg(vertices: &[P2], view_box: &str) -> std::io::Result<()> {
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
		"<polygon stroke-width='0.01%' stroke='black'  fill='none' points='"
	)?;
	for [x, y] in vertices {
		writeln!(file, "{},{}", x, y)?;
	}
	writeln!(file, "' />")?;
	for (i, [x, y]) in vertices.iter().enumerate() {
		let angle = rand::random::<f64>() * 6.28;
		let offset = 1.0;
		writeln!(
			file,
			"<text x='{}' y='{}' font-size='20%'>{}</text>",
			x + offset * angle.cos(),
			y + offset * angle.sin(),
			i
		)?;
	}
	// writeln!(
	// 	file,
	// 	"<circle r='1' cx='85.33333333333333' cy='42.666666666666664' fill='red'/>"
	// );
	writeln!(file, "</svg>")?;
	Ok(())
}
