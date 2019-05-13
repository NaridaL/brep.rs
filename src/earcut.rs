use bit_vec::BitVec;
use itertools::Itertools;
use num_traits::real::Real;
use std::fs::File;
use std::io::Write;
use std::iter::FromIterator;

pub fn earcut(data: &[f64], hole_indices: &[usize], dim: u8) {
	let has_holes = !hole_indices.is_empty();
	// let outerLen = holeIndices.get(0).unwrap_or(data.len()*dim);
}

fn earcut_no_holes(data: &[f64]) -> Vec<usize> {
	let modulo = |a: usize, b: usize| -> usize { ((a % b) + b) % b };
	let nv = data.len() / 2;
	let e = |i: usize| data[i];
	let mut triangles: Vec<usize> = Vec::with_capacity(3 * (nv - 2));
	let convex = |i, j, k| {
		triangle_is_convex(
			e(i * 2),
			e(i * 2 + 1),
			e(j * 2),
			e(j * 2 + 1),
			e(k * 2),
			e(k * 2 + 1),
		)
	};
	let mut is_convex = BitVec::from_fn(nv, |i| convex((i + nv - 1) % nv, i, (i + 1) % nv));
	println!("is_convex {:?}", is_convex);
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

	fn is_ear<F>(is_convex: &BitVec, is_used: &BitVec, e: F, i: usize, j: usize, k: usize) -> bool
	where
		F: Fn(usize) -> f64,
	{
		// println!("is_convex {:?}", is_convex[j]);
		is_convex[j]
			&& is_used.iter().enumerate().all(|(vi, used)| {
				used || vi == i
					|| vi == j || vi == k
					|| !triangle_contains_point(
						e(i * 2),
						e(i * 2 + 1),
						e(j * 2),
						e(j * 2 + 1),
						e(k * 2),
						e(k * 2 + 1),
						e(vi * 2),
						e(vi * 2 + 1),
					)
			})
	};

	fn is_ear_zorder<F, G>(
		is_convex: &BitVec,
		is_used: &BitVec,
		z_order_values: &[(usize, u32)],
		z_order_val: G,
		z_index_arr: &[usize],
		e: F,
		i: usize,
		j: usize,
		k: usize,
	) -> bool
	where
		F: Fn(usize) -> f64,
		G: Fn((f64, f64)) -> u32,
	{
		if !is_convex[j] {
			return false;
		}

		// println!("is_ear_zorder {} {} {}", i, j, k);

		let nv = z_order_values.len();

		let ax = e(i * 2);
		let ay = e(i * 2 + 1);
		let bx = e(j * 2);
		let by = e(j * 2 + 1);
		let cx = e(k * 2);
		let cy = e(k * 2 + 1);

		let min_triang_x = ax.min(bx).min(cx);
		let max_triang_x = ax.max(bx).max(cx);
		let min_triang_y = ay.min(by).min(cy);
		let max_triang_y = ay.max(by).max(cy);

		let min_triang_z = z_order_val((min_triang_x, min_triang_y));
		let max_triang_z = z_order_val((max_triang_x, max_triang_y));
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
		let result = a().interleave(b()).all(|(vi, _)| {
			// tested.push(vi);
			is_used[vi]
				|| vi == i || vi == j
				|| vi == k || !triangle_contains_point(ax, ay, bx, by, cx, cy, e(vi * 2), e(vi * 2 + 1))
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
							println!("{} {}", e(vi * 2), e(vi * 2 + 1));
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
			is_convex.set(i, is_convex[i] || convex(h, i, j));
			is_convex.set(j, is_convex[j] || convex(i, j, k));
		}
		current_indexes[0] = current_indexes[1];
		current_indexes[1] = current_indexes[2];
		current_indexes[2] = current_indexes[3];
		current_indexes[3] = unused_index_after(&next_unused_index, current_indexes[2]);
	}
	triangles
}

fn triangle_is_convex(ax: f64, ay: f64, bx: f64, by: f64, cx: f64, cy: f64) -> bool {
	// println!("{} {} {} {} {} {} ", ax, ay, bx, by, cx, cy);
	// println!("{}", (bx - ax) * (cy - ay) - (by - ay) * (cx - ax));
	(bx - ax) * (cy - ay) - (by - ay) * (cx - ax) > 0.0
}

fn triangle_contains_point(
	ax: f64,
	ay: f64,
	bx: f64,
	by: f64,
	cx: f64,
	cy: f64,
	px: f64,
	py: f64,
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

fn ddd<F>(f: F, v: P2) -> (f64, f64)
where
	F: Fn(f64, f64) -> f64,
{
	let eps = 1e-6;
	let (x, y) = v;
	let fxy = f(x, y);
	((f(x + eps, y) - fxy) / eps, (f(x, y + eps) - fxy) / eps)
}

fn curve_point<F>(f: &F, start: (f64, f64)) -> (f64, f64)
where
	F: Fn(f64, f64) -> f64,
{
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

fn interleave(a: u16, b: u16) -> u32 {
	let mut x: u32 = a.into();
	let mut y: u32 = b.into();

	x = (x | (x << 8)) & 0x00FF00FF;
	x = (x | (x << 4)) & 0x0F0F0F0F;
	x = (x | (x << 2)) & 0x33333333;
	x = (x | (x << 1)) & 0x55555555;

	y = (y | (y << 8)) & 0x00FF00FF;
	y = (y | (y << 4)) & 0x0F0F0F0F;
	y = (y | (y << 2)) & 0x33333333;
	y = (y | (y << 1)) & 0x55555555;

	x | (y << 1)
}

fn followAlgorithm2d<F>(
	ic: &F,
	start_p: (f64, f64),
	step_length: f64,
) -> (Vec<(f64, f64)>, Vec<(f64, f64)>)
where
	F: Fn(f64, f64) -> f64,
{
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

	#[test]
	fn tt3() {
		let f = |x: f64, y: f64| foo(x, -y, 3.0, 0.5);
		let start_point = curve_point(&f, (3.0, -1.0));
		let (points, _) = followAlgorithm2d(&f, start_point, 0.001);
		let (hole1_points, _) = followAlgorithm2d(&f, curve_point(&f, (0.1, 1.0)), 0.001);
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
		write_point_svg(&vertices, "-8 -8 16 16").unwrap();
		let triangles;
		{
			time_test!();
			triangles = earcut_no_holes(&vertices);
		}
		write_svg(&vertices, &triangles, "-8 -8 16 16").unwrap();
	}
}
fn lerp(a: f64, b: f64, t: f64) -> f64 {
	(1.0 - t) * a + t * b
}
use rand::Rng;
fn write_svg(vertices: &[f64], triangles: &[usize], view_box: &str) -> std::io::Result<()> {
	return Ok(());
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
	for (i, j, k) in triangles.iter().tuples() {
		writeln!(
			file,
			"<polygon points='{},{} {},{} {},{}' style='fill:rgb({},{},{});' />",
			vertices[i * 2],
			vertices[i * 2 + 1],
			vertices[j * 2],
			vertices[j * 2 + 1],
			vertices[k * 2],
			vertices[k * 2 + 1],
			rand::thread_rng().gen::<u8>(),
			rand::thread_rng().gen::<u8>(),
			rand::thread_rng().gen::<u8>(),
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
fn write_point_svg(vertices: &[f64], view_box: &str) -> std::io::Result<()> {
	return Ok(());
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
	for i in 0..vertices.len() / 2 {
		writeln!(file, "{},{} ", vertices[i * 2], vertices[i * 2 + 1])?;
	}
	writeln!(file, "' /></svg>")?;
	Ok(())
}
