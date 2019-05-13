#[macro_use]
extern crate time_test;

mod m4;
// mod NURBS;
mod earcut;
mod v4;

// #[cfg(test)]
// mod tests {
// 	// Note this useful idiom: importing names from outer (for mod tests) scope.
// 	use super::*;

// 	#[test]
// 	fn test_v4_add_assign() {
// 		let i: M4<i32> = M4::identity();
// 		let j = i + i;

// 		assert_eq!(j.0[0], 2);
// 	}

// 	// #[test]
// 	// fn test_bad_add() {
// 	//     // This assert would fire and test will fail.
// 	//     // Please note, that private functions can be tested too!
// 	//     assert_eq!(bad_add(1, 2), 3);
// 	// }
// }
